"""Command-line interface.

    clavagen run -q query.fa -r reference.fa -a annotation.json -o results/
    clavagen align -q query.fa -r reference.fa
    clavagen translate -s ATGGTG...
    clavagen example -o demo/        # write demo inputs and run them
    clavagen wizard                  # interactive prompt-driven input
    clavagen app                     # application window (what the start files open)
"""

from __future__ import annotations

import argparse
import json
import os
import sys

from . import __version__
from .align import ScoringScheme, align
from .genetic_code import find_orfs, gc_content, reverse_complement, translate
from .pipeline import Pipeline
from .report import write_all
from .seqio import load_sequences


def _add_run_args(p):
    p.add_argument("-q", "--query", action="append", required=True,
                   help="Query sequence(s): FASTA file, NAME=SEQUENCE, or raw sequence. Repeatable; multi-FASTA allowed.")
    p.add_argument("-r", "--reference", action="append", required=True,
                   help="Reference sequence(s) (DNA cDNA/genomic or protein). Repeatable; multi-FASTA allowed.")
    p.add_argument("-a", "--annotation", help="Annotation JSON (gene, CDS/exons, domains, known variants, population, in-silico).")
    p.add_argument("-g", "--gene", help="Gene symbol (uses built-in gene-disease knowledge if no annotation).")
    p.add_argument("--cds", help="CDS range on the reference, 1-based inclusive (e.g. 51-494). Overrides annotation.")
    p.add_argument("--pair", action="append", help="Explicit pairing QUERY_ID:REF_ID (otherwise best-scoring reference).")
    p.add_argument("--orthologs", action="append", help="Ortholog protein FASTA for conservation analysis.")
    p.add_argument("-e", "--evidence", action="append",
                   help="Reviewer evidence, e.g. 'c.20A>T=PS3,PP1_Strong', 'p.Glu7Val=-PM2', '*=PP4'. Repeatable.")
    p.add_argument("--evidence-file", help="JSON {variant: [criteria,...]} with reviewer evidence.")
    p.add_argument("-k", "--known", action="append",
                   help="Known pathogenic variants: a .txt/.csv/.tsv/.json file, or an entry like "
                        "'c.20A>T|p.Glu7Val|Pathogenic|source'. Repeatable. Added to the online ClinVar list.")
    p.add_argument("--offline", action="store_true",
                   help="Do not search online databases (ClinVar, RefSeq, UniProt, gnomAD, Ensembl VEP).")
    p.add_argument("--no-clinvar-class", action="store_true",
                   help="Report the automated ACMG class only (default: report the ClinVar class for known variants).")
    p.add_argument("--clinvar-min-stars", type=int, default=0, choices=[0, 1, 2, 3, 4],
                   help="Minimum ClinVar review stars for using its classification (default 0 = any record).")
    p.add_argument("--ref-setting", action="append",
                   help="Per-reference setting by reference ID or gene symbol, e.g. NM_000454.5:inheritance=AD/AR or SOD1:inheritance=AD/AR (repeatable).")
    p.add_argument("-p", "--phenotype", action="append", help="Patient phenotype term / HPO label (repeatable) for PP4.")
    p.add_argument("--clinical-evidence", metavar="JSON",
                   help="Clinical evidence you observed (segregation, phase, functional assays, de novo, case data...), "
                        "converted to ACMG criteria with ClinGen rules. See 'clavagen evidence-template'.")
    p.add_argument("--rules", metavar="JSON",
                   help="Gene rule profile (e.g. from a ClinGen expert panel specification): thresholds, criteria not "
                        "used, fixed strengths and caps. See 'clavagen evidence-template'.")
    p.add_argument("--patient-id"), p.add_argument("--sex", choices=["male", "female", "unknown"])
    p.add_argument("--age"), p.add_argument("--indication", help="Clinical indication for testing.")
    p.add_argument("-t", "--type", default="auto", choices=["auto", "dna", "protein"], help="Force query molecule type.")
    p.add_argument("--mode", default="local", choices=["local", "global"], help="Alignment mode (default local).")
    p.add_argument("--aligner", default="auto", choices=["auto", "python", "biopython"])
    p.add_argument("--strict", action="store_true", help="Disable heuristic (non-calibrated) evidence (in-silico heuristic, PP4 keyword match).")
    p.add_argument("--min-identity", type=float, default=80.0, help="QC warning threshold for alignment identity (%%).")
    p.add_argument("-o", "--outdir", default="clavagen_results")
    p.add_argument("--formats", default="tsv,vcf,json,md,html,svg,aln", help="Comma-separated output formats.")


def _evidence(args):
    ev = list(args.evidence or [])
    if getattr(args, "evidence_file", None):
        with open(args.evidence_file) as fh:
            for k, v in json.load(fh).items():
                ev.append(f"{k}={','.join(v) if isinstance(v, list) else v}")
    return ev


def _apply_rules(ann, path):
    """Merge a gene rule profile into the annotation. The file is either one profile
    {"source", "thresholds", "disabled", "strength", "max_strength", "start_loss_pvs1", "pp1_meioses"} for all
    references, or {"GENE or reference id": profile, ...}."""
    if not path:
        return ann
    with open(path) as fh:
        prof = json.load(fh)
    ann = dict(ann or {})
    keys = {"source", "thresholds", "disabled", "strength", "max_strength", "start_loss_pvs1", "pp1_meioses"}
    prof = {k: v for k, v in prof.items() if not k.startswith("_")}
    per = prof if not (set(prof) & keys) else {"*": prof}

    def put(target, p):
        p = {k: v for k, v in p.items() if not k.startswith("_")}
        if "thresholds" in p:
            target["thresholds"] = dict(target.get("thresholds") or {}, **p.pop("thresholds"))
        target["rules"] = dict(target.get("rules") or {}, **p)
    for k, p in per.items():
        if k == "*":
            put(ann, p)
        else:
            refs = ann.setdefault("references", {})
            put(refs.setdefault(k, {}), p)
    return ann


TEMPLATE = {
    "clinical_evidence": {
        "_help": "Keys: a variant (c. or p.) or '*' for all variants of this patient. Leave out what you do not have.",
        "c.20A>T": {
            "segregation": {"informative_meioses": 5, "families": 2},
            "in_trans": [{"partner": "Pathogenic", "phase": "confirmed"}],
            "homozygous_probands": 0,
            "functional": {"result": "damaging", "oddspath": 20, "assay": "enzyme activity <10%"},
            "de_novo": [{"confirmed": True, "phenotype": "specific"}],
            "case_control": {"odds_ratio": 6.2, "ci_lower": 2.1},
            "affected_probands": {"count": 4, "strength": "Moderate"},
            "phenotype_specific": True,
            "healthy_adults": {"count": 0, "zygosity": "homozygous"},
            "non_segregation": False, "benign_phase": False, "alternative_cause": False
        }
    },
    "rules": {
        "_help": "One profile for all references, or {'GENE': profile}. Copy values from the expert panel "
                 "specification (ClinGen CSpec Registry, https://cspec.clinicalgenome.org) - do not guess them.",
        "MYGENE": {
            "source": "ClinGen <panel> VCEP specification v<version> (<CSpec URL>)",
            "thresholds": {"ba1": 0.001, "bs1": 0.0002, "pm2": 0.00002},
            "disabled": ["PP2"],
            "strength": {"PM2": "Supporting"},
            "max_strength": {"PP3": "Supporting"},
            "start_loss_pvs1": "Moderate",
            "pp1_meioses": {"Supporting": 3, "Moderate": 5, "Strong": 7}
        }
    }
}


def cmd_evidence_template(args):
    out = args.output or "clavagen_evidence_template.json"
    with open(out, "w") as fh:
        json.dump(TEMPLATE, fh, indent=2)
    print(f"Template written to {os.path.abspath(out)}. Save the 'clinical_evidence' part as one file "
          "(--clinical-evidence) and the 'rules' part as another (--rules), and replace the example values.")


def cmd_run(args):
    patient = {"id": args.patient_id, "sex": args.sex, "age": args.age, "indication": args.indication}
    patient = {k: v for k, v in patient.items() if v}
    from .online import parse_known_variants_text
    known = []
    for k in (getattr(args, "known", None) or []):
        if os.path.isfile(k):
            with open(k) as fh:
                known += parse_known_variants_text(fh.read(), default_source=os.path.basename(k))
        else:
            known += parse_known_variants_text(k, default_source="Command line")
    online = not getattr(args, "offline", False)
    from .annotation import load_annotation
    ann = _apply_ref_settings(load_annotation(args.annotation) if args.annotation else None, getattr(args, "ref_setting", None))
    ann = _apply_rules(ann, getattr(args, "rules", None))
    clin = None
    if getattr(args, "clinical_evidence", None):
        with open(args.clinical_evidence) as fh:
            clin = json.load(fh)
    pl = Pipeline(args.query, args.reference, annotation=ann, gene=args.gene, cds=args.cds,
                  online=online, known_variants=known, log=lambda m: print("  ..", m),
                  use_clinvar=not getattr(args, "no_clinvar_class", False),
                  clinvar_min_stars=getattr(args, "clinvar_min_stars", 0),
                  pairs=args.pair, orthologs=args.orthologs, evidence=_evidence(args), phenotypes=args.phenotype,
                  mol_type=args.type, aligner=args.aligner, mode=args.mode, strict=args.strict,
                  patient=patient, min_identity=args.min_identity, clinical_evidence=clin)
    pl.annotation_path = args.annotation
    res = pl.run()
    files = write_all(res, args.outdir, tuple(args.formats.split(",")))
    _print_summary(res)
    print(f"\nOutputs written to {os.path.abspath(args.outdir)}:")
    for f in files:
        print("  ", os.path.basename(f))
    return res


def _print_summary(res):
    for s in res["samples"]:
        st = s["alignment"].stats
        print(f"\n== {s['query'].id} vs {s['reference'].id}  [{s['mode']}]  identity {st['identity_pct']}%  "
              f"gene {s['context'].get('gene', '-')}")
        for w in s["qc"].get("warnings", []):
            print(f"   ! {w}")
        if not s["variants"]:
            print("   no variants")
        for v in s["variants"]:
            print(f"   {v['id']:7} {v['hgvs_c_short'] or '':24} {v['hgvs_p'] or '':26} {v['consequence_label']:18} "
                  f"{v['zygosity_short']:14} {v['classification']:24} [{', '.join(v['criteria_met'])}]")
        print(f"   => {s['clinical']['status']}")


def _apply_ref_settings(annotation, settings):
    """--ref-setting NM_000454.5:inheritance=AD/AR or SOD1:inheritance=AD/AR - by reference ID or gene
    symbol (also gene=, cds=51-494, transcript=)."""
    if not settings:
        return annotation
    ann = dict(annotation or {})
    refs = dict(ann.get("references") or {})
    for item in settings:
        rid, _, kv = item.partition(":")
        key, _, val = kv.partition("=")
        if not rid or not key:
            raise ValueError(f"--ref-setting must look like REF_ID:key=value (got '{item}')")
        val = [int(x) for x in val.replace("..", "-").split("-")] if key.strip() == "cds" else val.strip()
        refs.setdefault(rid.strip(), {})[key.strip()] = val
    ann["references"] = refs
    return ann


def cmd_cohort(args):
    from .annotation import load_annotation
    from .cohort import read_patient_table, run_cohort
    from .online import parse_known_variants_text
    patients = read_patient_table(args.patients)
    known = []
    for k in args.known or []:
        if os.path.isfile(k):
            with open(k) as fh:
                known += parse_known_variants_text(fh.read(), default_source=os.path.basename(k))
        else:
            known += parse_known_variants_text(k, default_source="Command line")
    ref = args.reference[0] if len(args.reference) == 1 else "\n".join(args.reference)
    if len(args.reference) == 1 and os.path.isfile(args.reference[0]):
        ref = args.reference[0]
    ann = _apply_ref_settings(load_annotation(args.annotation) if args.annotation else None, getattr(args, "ref_setting", None))
    ref = list(args.reference)
    summ = run_cohort(patients, ref, args.outdir, annotation=ann,
                      gene=args.gene, cds=args.cds, strict=args.strict, online=not args.offline, known=known,
                      use_clinvar=not args.no_clinvar_class, clinvar_min_stars=args.clinvar_min_stars,
                      log=lambda m: print("  ", m))
    folder = os.path.join(os.path.abspath(args.outdir), summ["name"])
    print(f"\n{len(summ['patients'])} patient(s) analysed. Results: {folder}")
    print("  cohort_summary.html / cohort_summary.tsv + one sub-folder per patient")
    for e in summ["patients"]:
        print(f"  {'OK ' if e['ok'] else 'ERR'} {e['id']}" + ("" if e["ok"] else f": {e.get('error')}"))


def cmd_benchmark(args):
    import datetime as _dt
    import hashlib
    import platform
    import sys
    from .benchmark import EREPO_URL, load_erepo, run_benchmark
    meta = {"command": " ".join(sys.argv), "python": platform.python_version(),
            "platform": f"{platform.system()} {platform.release()} ({platform.machine()})"}
    if args.erepo:
        with open(args.erepo, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        meta["erepo_file"] = os.path.basename(args.erepo)
        meta["erepo_date"] = _dt.datetime.fromtimestamp(os.path.getmtime(args.erepo)).date().isoformat()
    else:
        from .online import http
        print("Downloading the ClinGen Evidence Repository export...")
        text = http(EREPO_URL, kind="text", timeout=120)
        os.makedirs(args.outdir, exist_ok=True)
        with open(os.path.join(args.outdir, "erepo_download.tsv"), "w") as fh:
            fh.write(text)
        meta["erepo_file"] = "erepo_download.tsv"
        meta["erepo_date"] = _dt.date.today().isoformat()
    meta["erepo_sha256"] = hashlib.sha256(text.encode("utf-8")).hexdigest()
    meta["erepo_url"] = EREPO_URL
    try:
        from .online import http
        meta["ensembl_release"] = (http("https://rest.ensembl.org/info/software?content-type=application/json",
                                        use_cache=False, timeout=15) or {}).get("release")
    except Exception:
        meta["ensembl_release"] = None
    recs, stats = load_erepo(text)
    meta["erepo_filter"] = stats
    print(f"Evidence Repository: {stats['total']} rows; usable exonic NM_ variants: {stats['used']} "
          f"(retracted {stats['retracted']}, no transcript HGVS {stats['no_transcript_hgvs']}, "
          f"intronic/complex {stats['intronic_or_complex']})")
    genes = [g.strip() for g in args.genes.split(",")] if args.genes else None
    summ = run_benchmark(recs, args.outdir, flank=args.flank, max_n=args.max or None, seed=args.seed, genes=genes,
                         resume=args.resume, meta=meta, figures=not args.no_figures)
    pct = lambda x: "-" if x is None else f"{100 * x:.1f}%"  # noqa: E731
    print(f"\nAnalysed {summ['analysed']} / {summ['attempted']} variants")
    print(f"  3-tier agreement (P/LP, VUS, B/LB): {pct(summ['three_tier_agreement'])}")
    print(f"  exact 5-class agreement:            {pct(summ['exact_agreement'])}")
    print(f"  critical discordance (P/LP<->B/LB): {pct(summ['critical_discordance'])}")
    k = (summ.get("statistics") or {}).get("kappa_three_tier") or {}
    if k.get("kappa") is not None:
        print(f"  Cohen's kappa (3-tier):             {k['kappa']:.3f} (95% CI {k['ci95'][0]}-{k['ci95'][1]})")
    print(f"Report: {os.path.join(os.path.abspath(args.outdir), 'benchmark_report.html')}")
    print(f"Tables, methods text and figures: {os.path.abspath(args.outdir)}")


def cmd_benchmark_figures(args):
    import json
    from .benchmark import summarise
    from .benchmark_figures import make_figures
    d = args.folder
    with open(os.path.join(d, "benchmark_summary.json")) as fh:
        old = json.load(fh)
    with open(os.path.join(d, "benchmark_results.jsonl")) as fh:
        results = [json.loads(line) for line in fh if line.strip()]
    summ = summarise(results, meta=old.get("run"), seed=(old.get("run") or {}).get("seed", 7))
    from .benchmark import write_outputs
    write_outputs(results, summ, d, figures=False)      # summary, tables, methods text and report, recounted
    made = make_figures(summ, results, os.path.join(d, "figures"))
    print(f"Summary, tables and report updated; {len(made)} figure files written to "
          f"{os.path.join(os.path.abspath(d), 'figures')}")


def cmd_check(args):
    from .online import check_connection
    for name, ok, msg in check_connection():
        print(f"  {'OK  ' if ok else 'FAIL'} {name}: {msg}")


def cmd_align(args):
    qs, rs = load_sequences(args.query, "query"), load_sequences(args.reference, "ref")
    for q in qs:
        for r in rs:
            mol = "protein" if "protein" in (q.mol_type, r.mol_type) else "dna"
            a = align(r.id, r.seq, q.id, q.seq, mol_type=mol, mode=args.mode, backend=args.aligner)
            print(a.format_text(ScoringScheme.default(mol)))


def cmd_translate(args):
    recs = load_sequences(args.sequence, "seq", "dna")
    for r in recs:
        s = reverse_complement(r.seq) if args.reverse else r.seq
        print(f">{r.id} length={len(s)} GC={gc_content(s):.1f}%")
        for f in range(3):
            print(f"frame +{f + 1}: {translate(s[f:])}")
        orfs = find_orfs(s, args.min_orf)
        for o in orfs[:5]:
            print(f"ORF {o[0] + 1}-{o[1]} ({len(o[2]) - 1} aa): {o[2]}")


def cmd_example(args):
    from .examples import write_examples
    out = args.outdir
    paths = write_examples(os.path.join(out, "inputs"))
    print("Demo inputs written:", *paths.values(), sep="\n  ")
    ns = argparse.Namespace(
        query=[paths["dna_queries"], paths["protein_query"]], reference=[paths["reference_dna"]],
        annotation=paths["annotation"], gene=None, cds=None, pair=None, orthologs=None,
        evidence=["c.20A>T=PS3,PS4"], evidence_file=None,
        phenotype=["microcytic anaemia", "vaso-occlusive pain crises"],
        patient_id="DEMO-001", sex="female", age="24", indication="Anaemia work-up (demonstration data)",
        type="auto", aligner="auto", mode="local", strict=False, min_identity=80.0, known=None, offline=True,
        outdir=os.path.join(out, "results"), formats="tsv,vcf,json,md,html,svg,aln")
    return cmd_run(ns)


def cmd_wizard(args):
    print("CLAVAGEN interactive input. Paste sequences (FASTA or raw); finish each with an empty line.\n")

    def read_block(prompt):
        print(prompt)
        lines = []
        while True:
            try:
                ln = input()
            except EOFError:
                break
            if not ln.strip():
                break
            lines.append(ln)
        return "\n".join(lines)

    def as_spec(text, role):
        if not text.strip():
            return []
        if os.path.isfile(text.strip()):
            return [text.strip()]
        import tempfile
        fd, p = tempfile.mkstemp(suffix=".fa", prefix=role)
        with os.fdopen(fd, "w") as fh:
            fh.write(text if text.lstrip().startswith(">") else f">{role}1\n{text}\n")
        return [p]

    q = as_spec(read_block("Query sequence(s) or file path:"), "query")
    r = as_spec(read_block("Reference sequence(s) or file path:"), "reference")
    ann = input("Annotation JSON path (optional): ").strip() or None
    gene = input("Gene symbol (optional): ").strip() or None
    ph = [x.strip() for x in input("Phenotypes, comma-separated (optional): ").split(",") if x.strip()]
    ev = [x.strip() for x in input("Reviewer evidence e.g. c.20A>T=PS3 (optional, ';'-separated): ").split(";") if x.strip()]
    out = input("Output folder [clavagen_results]: ").strip() or "clavagen_results"
    ns = argparse.Namespace(query=q, reference=r, annotation=ann, gene=gene, cds=None, pair=None, orthologs=None,
                            evidence=ev, evidence_file=None, phenotype=ph, patient_id=None, sex=None, age=None,
                            indication=None, type="auto", aligner="auto", mode="local", strict=False,
                            min_identity=80.0, outdir=out, formats="tsv,vcf,json,md,html,svg,aln",
                            known=None, offline=False)
    cmd_run(ns)


def cmd_app(args):
    from .webui import main as web_main
    web_main(out_base=getattr(args, "outdir", None), port=getattr(args, "port", None),
             open_browser=True, window=getattr(args, "window", "auto"))


def cmd_web(args):
    from .webui import main as web_main
    web_main(out_base=args.outdir, port=args.port, open_browser=not args.no_browser, window="browser")


def cmd_tk(args):
    from .gui import main as gui_main
    gui_main()


def main(argv=None):
    ap = argparse.ArgumentParser(prog="clavagen", description="Clinical variant calling and ACMG/AMP classification "
                                 "from DNA or protein sequences.")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="Call, annotate and classify variants; write maps, alignments and reports.")
    _add_run_args(p)
    p.set_defaults(func=cmd_run)
    p = sub.add_parser("cohort", help="Analyse many patients from a patient table (CSV/TSV/JSON).")
    p.add_argument("--patients", required=True,
                   help="Table with columns id, sex, age, indication, phenotypes (comma-sep), evidence (;-sep), "
                        "query (sequence or FASTA file path).")
    p.add_argument("-r", "--reference", action="append", required=True, help="Reference file, sequence or NM_ accession.")
    p.add_argument("-a", "--annotation"), p.add_argument("-g", "--gene"), p.add_argument("--cds")
    p.add_argument("-k", "--known", action="append", help="Known pathogenic variants file or entry (repeatable).")
    p.add_argument("--offline", action="store_true"), p.add_argument("--strict", action="store_true")
    p.add_argument("--no-clinvar-class", action="store_true", help="Report the automated ACMG class only.")
    p.add_argument("--clinvar-min-stars", type=int, default=0, choices=[0, 1, 2, 3, 4])
    p.add_argument("--ref-setting", action="append",
                   help="Per-reference setting, e.g. NM_000454.5:inheritance=AD/AR, NM_000518.5:cds=51-494, REF:gene=HBB")
    p.add_argument("-o", "--outdir", default="clavagen_results")
    p.set_defaults(func=cmd_cohort)
    p = sub.add_parser("benchmark", help="Benchmark the ACMG engine against ClinGen Evidence Repository expert classifications.")
    p.add_argument("--erepo", help="Evidence Repository tabbed export (TSV). Omit to download it (internet needed).")
    p.add_argument("--max", type=int, default=200, help="Maximum variants (stratified sample; default 200, 0 = all).")
    p.add_argument("--genes", help="Comma-separated gene list to restrict the benchmark.")
    p.add_argument("--flank", type=int, default=150, help="Bases either side of the variant in the simulated read.")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--resume", action="store_true",
                   help="Continue an interrupted run in the same output folder (same --erepo, --max, --genes, --seed).")
    p.add_argument("--no-figures", action="store_true", help="Skip the publication figures.")
    p.add_argument("-o", "--outdir", default="clavagen_benchmark")
    p.set_defaults(func=cmd_benchmark)
    p = sub.add_parser("evidence-template", help="Write an example clinical-evidence and gene-rules JSON file.")
    p.add_argument("-o", "--output")
    p.set_defaults(func=cmd_evidence_template)
    p = sub.add_parser("benchmark-figures", help="Recount the summary/tables and re-draw the figures of a finished benchmark run (no re-analysis).")
    p.add_argument("folder", help="Benchmark output folder (with benchmark_summary.json and benchmark_results.jsonl).")
    p.set_defaults(func=cmd_benchmark_figures)
    p = sub.add_parser("check-connection", help="Test access to NCBI, UniProt, gnomAD and Ensembl.")
    p.set_defaults(func=cmd_check)
    p = sub.add_parser("align", help="Local/global alignment only.")
    p.add_argument("-q", "--query", action="append", required=True)
    p.add_argument("-r", "--reference", action="append", required=True)
    p.add_argument("--mode", default="local", choices=["local", "global"])
    p.add_argument("--aligner", default="auto", choices=["auto", "python", "biopython"])
    p.set_defaults(func=cmd_align)
    p = sub.add_parser("translate", help="Six/three-frame translation, ORFs and GC content.")
    p.add_argument("-s", "--sequence", action="append", required=True)
    p.add_argument("--reverse", action="store_true", help="Use reverse complement.")
    p.add_argument("--min-orf", type=int, default=30)
    p.set_defaults(func=cmd_translate)
    p = sub.add_parser("example", help="Write a worked HBB demonstration dataset and run it.")
    p.add_argument("-o", "--outdir", default="clavagen_demo")
    p.set_defaults(func=cmd_example)
    p = sub.add_parser("wizard", help="Interactive prompt-driven input.")
    p.set_defaults(func=cmd_wizard)
    for name in ("app", "gui"):
        p = sub.add_parser(name, help="Open the CLAVAGEN application window." if name == "app" else "Alias of 'app'.")
        p.add_argument("-o", "--outdir", help="Folder for results (default ./results).")
        p.add_argument("--port", type=int)
        p.add_argument("--window", default="auto", choices=["auto", "native", "app", "browser"],
                       help="native = pywebview window; app = Edge/Chrome app window; browser = normal tab.")
        p.set_defaults(func=cmd_app)
    p = sub.add_parser("web", help="Serve the interface in a normal browser tab.")
    p.add_argument("-o", "--outdir", help="Folder for results (default ./results).")
    p.add_argument("--port", type=int)
    p.add_argument("--no-browser", action="store_true", help="Only start the local server.")
    p.set_defaults(func=cmd_web)
    p = sub.add_parser("tk", help="Minimal tkinter window (legacy).")
    p.set_defaults(func=cmd_tk)
    args = ap.parse_args(argv)
    try:
        args.func(args)
    except (ValueError, FileNotFoundError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
