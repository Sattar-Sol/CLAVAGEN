"""Multi-patient (cohort) analysis.

Each patient has their own query sequence(s), demographics, phenotypes and reviewer
evidence; the reference, annotation, known variants and settings are shared. Output:

<run folder>/
    run_summary.json             overview used by the application window / history
    cohort_summary.tsv           one row per variant, all patients
    cohort_summary.html          overview table of all patients
    <patient_id>/                full per-patient outputs (report.html, maps, alignments, ...)
        inputs/
"""

from __future__ import annotations

import csv
import datetime as _dt
import io
import json
import os
import re
from html import escape

from . import __version__

CLASS_ORDER = ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]
CLASS_COLORS = {"Pathogenic": "#b2182b", "Likely pathogenic": "#ef6548", "Uncertain significance": "#e3a008",
                "Likely benign": "#6baed6", "Benign": "#2166ac"}
ACCESSION = re.compile(r"N[MR]_\d+(\.\d+)?")


def safe(s):
    return re.sub(r"[^\w.\-]+", "_", s).strip("_") or "patient"


def _split_list(v, sep):
    if isinstance(v, (list, tuple)):
        return [str(x).strip() for x in v if str(x).strip()]
    return [x.strip() for x in str(v or "").split(sep) if x.strip()]


def normalise_patients(patients):
    """Validate and normalise patient dicts. Required: id, query (sequence text or file path)."""
    out, seen = [], set()
    for i, p in enumerate(patients or [], 1):
        pid = str(p.get("id") or "").strip() or f"Patient{i}"
        if pid in seen:
            raise ValueError(f"Duplicate patient ID: {pid}")
        seen.add(pid)
        qv = p.get("query") or ""
        q = [str(x).strip() for x in qv if str(x).strip()] if isinstance(qv, list) else str(qv).strip()
        traces = [t for t in (p.get("traces") or []) if t.get("data")]
        if not q and not traces:
            raise ValueError(f"Patient {pid}: no sequence or trace file given.")
        rel = str(p.get("relationship") or "").strip().lower()
        aff = str(p.get("affected") or "").strip().lower()
        aff = {"yes": "affected", "y": "affected", "1": "affected", "true": "affected", "no": "unaffected",
               "n": "unaffected", "0": "unaffected", "false": "unaffected"}.get(aff, aff)
        pc = p.get("parentage_confirmed")
        pc = pc if isinstance(pc, bool) else str(pc or "").strip().lower() in ("1", "yes", "y", "true")
        out.append({"id": pid, "sex": str(p.get("sex") or "").strip(), "age": str(p.get("age") or "").strip(),
                    "family": str(p.get("family") or "").strip(),
                    "reference": str(p.get("reference") or "").strip(),
                    "relationship": rel if rel in ("proband", "mother", "father", "sibling", "child", "partner", "other") else (
                        "other" if rel else ""),
                    "affected": aff if aff in ("affected", "unaffected") else "unknown",
                    "parentage_confirmed": pc,
                    "indication": str(p.get("indication") or "").strip(),
                    "phenotypes": _split_list(p.get("phenotypes"), ","),
                    "evidence": _split_list(p.get("evidence"), ";"),
                    "clinical_evidence": _clinical(p.get("clinical_evidence")),
                    "query": q, "traces": traces, "notes": str(p.get("notes") or "").strip()})
    if not out:
        raise ValueError("Add at least one patient.")
    return out


def _clinical(x):
    """Guided clinical evidence: a dict, a JSON string, or empty."""
    if not x:
        return {}
    if isinstance(x, dict):
        return x
    try:
        v = json.loads(x)
        return v if isinstance(v, dict) else {}
    except (TypeError, ValueError):
        raise ValueError(f"clinical_evidence is not valid JSON: {str(x)[:80]}")


def read_patient_table(path_or_text):
    """Read a patient list from CSV/TSV (header: id, sex, age, indication, phenotypes, evidence, query)
    or JSON. ``query`` may be a sequence or a path to a FASTA file (relative to the table)."""
    base = ""
    text = path_or_text
    if os.path.isfile(path_or_text):
        base = os.path.dirname(os.path.abspath(path_or_text))
        with open(path_or_text, encoding="utf-8-sig") as fh:
            text = fh.read()
    text = text.strip()
    if text.startswith("["):
        rows = json.loads(text)
    else:
        dialect = "excel-tab" if text.split("\n", 1)[0].count("\t") > text.split("\n", 1)[0].count(",") else "excel"
        rows = list(csv.DictReader(io.StringIO(text), dialect=dialect))
    pts = []
    for r in rows:
        r = {str(k).strip().lower(): (v or "").strip() if isinstance(v, str) else v for k, v in r.items() if k}
        q = r.get("query") or r.get("sequence") or r.get("fasta") or r.get("trace") or ""
        parts = [x.strip() for x in q.split(";")] if (";" in q and len(q) < 2000) else [q]
        resolved = []
        for part in parts:
            cand = os.path.join(base, part) if base and part and not os.path.isabs(part) else part
            if part and len(part) < 400 and os.path.isfile(cand):
                if open(cand, "rb").read(4) == b"ABIF":
                    resolved.append(os.path.abspath(cand))   # Sanger trace file
                else:
                    with open(cand) as fh:
                        resolved.append(fh.read())
            elif part:
                resolved.append(part)
        q = resolved if len(resolved) > 1 or (resolved and os.path.isfile(resolved[0])) else (resolved[0] if resolved else "")
        pts.append({"id": r.get("id") or r.get("patient_id") or r.get("patient"), "sex": r.get("sex", ""),
                    "reference": r.get("reference") or r.get("gene") or "",
                    "family": r.get("family") or r.get("family_id") or "", "relationship": r.get("relationship") or
                    r.get("relation") or "", "affected": r.get("affected", ""),
                    "parentage_confirmed": r.get("parentage_confirmed", ""),
                    "age": r.get("age", ""), "indication": r.get("indication", ""),
                    "phenotypes": r.get("phenotypes", ""), "evidence": r.get("evidence", ""),
                    "clinical_evidence": r.get("clinical_evidence", ""), "query": q})
    return normalise_patients(pts)


def _write_fasta_text(text, path, role):
    text = text.strip()
    if not text.startswith(">"):
        text = f">{role}\n{text}"
    with open(path, "w") as fh:
        fh.write(text + "\n")
    return path


def run_cohort(patients, reference, out_base, annotation=None, gene=None, cds=None, strict=False,
               online=True, known=None, run_name=None, log=None, formats=None, use_clinvar=True,
               clinvar_min_stars=0):
    """Analyse every patient against the shared reference(s). Returns the run summary dict."""
    from .pipeline import Pipeline
    from .report import write_all
    log = log or (lambda *_: None)
    patients = normalise_patients(patients)
    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    if not run_name:
        run_name = f"{stamp}_{safe(patients[0]['id'])}" if len(patients) == 1 else f"{stamp}_{len(patients)}_patients"
    folder = os.path.join(out_base, run_name)
    shared = os.path.join(folder, "shared_inputs")
    os.makedirs(shared, exist_ok=True)
    # shared inputs
    if isinstance(reference, str):
        ref_lines = [x.strip() for x in reference.strip().splitlines() if x.strip()]
        if ref_lines and all(ACCESSION.fullmatch(x) for x in ref_lines):
            refs = ref_lines
        elif os.path.isfile(reference.strip()):
            refs = [reference.strip()]
        else:
            refs = [_write_fasta_text(reference, os.path.join(shared, "reference.fa"), "reference")]
    else:
        # list of references: file paths, RefSeq accessions, or {"fasta": text} / {"accession": "NM_..."}
        refs = []
        for k, item in enumerate(reference, 1):
            if isinstance(item, dict):
                if item.get("accession"):
                    refs.append(item["accession"].strip())
                elif item.get("fasta"):
                    refs.append(_write_fasta_text(item["fasta"], os.path.join(shared, f"reference_{k}.fa"),
                                                  item.get("id") or f"reference{k}"))
            else:
                refs.append(item)
    if annotation:
        with open(os.path.join(shared, "annotation.json"), "w") as fh:
            json.dump(annotation, fh, indent=2)
    if known:
        with open(os.path.join(shared, "known_variants.json"), "w") as fh:
            json.dump(known, fh, indent=2)
    auto_cache = {}
    entries = []
    results, prepared = {}, {}

    def analyse(p, family_ev=None):
        patient_meta = {k: p[k] for k in ("id", "sex", "age", "indication") if p.get(k)}
        pl = Pipeline(prepared[p["id"]], refs, annotation=annotation, gene=gene, cds=cds, evidence=p["evidence"],
                      phenotypes=p["phenotypes"], strict=strict, patient=patient_meta, online=online,
                      known_variants=known, log=log, auto_cache=auto_cache, use_clinvar=use_clinvar,
                      clinvar_min_stars=clinvar_min_stars, family_evidence=family_ev,
                      default_reference=p.get("reference") or None,
                      clinical_evidence=p.get("clinical_evidence") or None)
        pl.annotation_path = "shared_inputs/annotation.json" if annotation else None
        return pl.run()

    for i, p in enumerate(patients, 1):
        log(f"[{i}/{len(patients)}] Analysing {p['id']}...")
        pdir = os.path.join(folder, safe(p["id"]))
        os.makedirs(os.path.join(pdir, "inputs"), exist_ok=True)
        specs = []
        qs = p["query"] if isinstance(p["query"], list) else ([p["query"]] if p["query"] else [])
        for k, q in enumerate(qs, 1):
            if len(q) < 400 and os.path.isfile(q):
                if open(q, "rb").read(4) == b"ABIF":
                    import shutil
                    dst = os.path.join(pdir, "inputs", os.path.basename(q))
                    shutil.copyfile(q, dst)
                    specs.append(dst)
                else:
                    specs.append(q)
            else:
                specs.append(_write_fasta_text(q, os.path.join(pdir, "inputs", f"query{'' if k == 1 else k}.fa"),
                                               safe(p["id"])))
        import base64
        for t in p.get("traces", []):
            fn = safe(os.path.splitext(t.get("name") or "trace")[0]) + ".ab1"
            dst = os.path.join(pdir, "inputs", fn)
            raw = t["data"].split(",", 1)[1] if "," in t["data"] else t["data"]
            with open(dst, "wb") as fh:
                fh.write(base64.b64decode(raw))
            specs.append(dst)
        prepared[p["id"]] = specs
        with open(os.path.join(pdir, "inputs", "patient.json"), "w") as fh:
            json.dump({k: v for k, v in p.items() if k not in ("query", "traces")}, fh, indent=2)
        try:
            results[p["id"]] = analyse(p)
        except Exception as e:  # one failing patient must not stop the others
            with open(os.path.join(pdir, "ERROR.txt"), "w") as fh:
                fh.write(str(e))
            entries.append({"id": p["id"], "folder": safe(p["id"]), "ok": False, "error": str(e)})
            log(f"   ! {p['id']}: {e}")

    # ---- family / couple analysis, then re-classify members with family evidence
    from .family import analyse_family, group_families, vkey
    fams = group_families(patients)
    fam_results = {}
    for fid, members in fams.items():
        log(f"Family analysis: {fid} ({len(members)} members)")
        fa = analyse_family(fid, members, results)
        if fa["evidence"]:
            for m in members:
                if m["id"] not in results:
                    continue
                keys = {vkey(v) for smp in results[m["id"]]["samples"] for v in smp["variants"]}
                if keys & set(fa["evidence"]):
                    try:
                        results[m["id"]] = analyse(m, family_ev=fa["evidence"])
                    except Exception as e:
                        log(f"   ! re-analysis of {m['id']} failed: {e}")
            fa = dict(analyse_family(fid, members, results), evidence=fa["evidence"])
        fam_results[fid] = fa
        for m in members:
            if m["id"] in results:
                results[m["id"]]["family"] = fa

    for p in patients:
        if p["id"] not in results:
            continue
        pdir = os.path.join(folder, safe(p["id"]))
        res = results[p["id"]]
        if res.get("family"):
            _family_to_clinical(res, p)
        write_all(res, pdir, formats or ("tsv", "vcf", "json", "md", "html", "svg", "aln"))
        entries.append({"id": p["id"], "folder": safe(p["id"]), "ok": True,
                        "family": p.get("family") or "", "relationship": p.get("relationship") or "",
                        "affected": p.get("affected") or ""})
    order = {p["id"]: i for i, p in enumerate(patients)}
    entries.sort(key=lambda e: order.get(e["id"], 0))
    if fam_results:
        from .family import family_html, pedigree_svg
        for fid, fa in fam_results.items():
            with open(os.path.join(folder, f"family_{safe(fid)}_pedigree.svg"), "w") as fh:
                fh.write(pedigree_svg(fa))
            from .report import svg_to_png
            svg_to_png(os.path.join(folder, f"family_{safe(fid)}_pedigree.svg"))
    summary = {"software": "CLAVAGEN", "version": __version__, "name": run_name,
               "date": _dt.datetime.now().isoformat(timespec="seconds"), "online": online,
               "patients": entries,
               "families": [_family_json(fa) for fa in fam_results.values()]}
    with open(os.path.join(folder, "run_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    write_cohort_tables(folder, summary)
    return summary


def _family_json(fa):
    d = {k: v for k, v in fa.items() if k != "evidence"}
    d["evidence"] = [{"key": k, "code": c, "strength": st, "rationale": why}
                     for k, lst in fa["evidence"].items() for c, st, why in lst]
    d["pedigree"] = f"family_{safe(fa['family'])}_pedigree.svg"
    return d


def _family_to_clinical(res, p):
    """Add family findings to the patient's clinical summary."""
    fa = res["family"]
    if fa.get("phases"):
        from .genotype import apply_patient_genotype
        apply_patient_genotype(res["samples"], {k: p.get(k) for k in ("id", "sex", "age", "indication")},
                               fa["phases"])
    for s in res["samples"]:
        c = s["clinical"]
        mine = [x for x in fa["findings"]]
        c["family"] = {"family": fa["family"], "relationship": p.get("relationship") or "",
                       "affected": p.get("affected") or "", "proband": fa["proband"], "findings": mine,
                       "couples": [cp for cp in fa["couples"] if p["id"] in cp["partners"]]}
        for cp in c["family"]["couples"]:
            for r in cp["risks"]:
                if r["level"] == "high":
                    c["recommendations"].insert(0, f"Reproductive risk ({cp['label']}): {r['risk']} - genetic "
                                                   "counselling and reproductive options.")
                    break


def patient_summary(folder, sub):
    """Summary of one patient's outputs (reads clinical_summary.json)."""
    pdir = os.path.join(folder, sub)
    with open(os.path.join(pdir, "clinical_summary.json")) as fh:
        js = json.load(fh)
    samples, allc = [], []
    for s in js["samples"]:
        b = safe_base(s)
        vs = s["variants"]
        allc += [v["classification"] for v in vs]
        samples.append({
            "query": s["query"]["id"], "reference": s["reference"]["id"],
            "status": s["clinical_summary"]["status"], "highest": s["clinical_summary"]["highest_class"],
            "warnings": s["qc"].get("warnings", []), "map": b + "_map.svg", "aln": b + "_alignment.txt",
            "variants": [{"id": v["id"], "c": v.get("hgvs_c_short") or "-", "p": v.get("hgvs_p") or "-",
                          "hgvs_c_full": v.get("hgvs_c_full") or "", "gene": v.get("gene") or "",
                          "type": v["consequence_label"], "zyg": v.get("zygosity_short", ""),
                          "cls": v["classification"], "pts": v.get("points"),
                          "acmg": v.get("acmg_classification", v["classification"]),
                          "basis": v.get("classification_basis", ""),
                          "clinvar": v.get("clinvar_classification", ""), "stars": v.get("clinvar_stars"),
                          "crit": ", ".join(v["criteria_met"]) or "-"} for v in vs]})
    pt = js["meta"].get("patient") or {}
    return {"patient": pt, "phenotypes": js["meta"].get("phenotypes", []), "samples": samples,
            "n_variants": len(allc), "highest": next((c for c in CLASS_ORDER if c in allc), None),
            "date": js["meta"]["run_date"]}


def safe_base(s):
    return re.sub(r"[^\w.\-]+", "_", f"{s['query']['id']}__vs__{s['reference']['id']}")


def write_cohort_tables(folder, summary):
    rows, prow = [], []
    for e in summary["patients"]:
        if not e["ok"]:
            prow.append((e, None))
            continue
        ps = patient_summary(folder, e["folder"])
        prow.append((e, ps))
        for s in ps["samples"]:
            if not s["variants"]:
                rows.append([e["id"], s["query"], "", "", "", "", "", "No variants detected", "", "", s["status"]])
            for v in s["variants"]:
                rows.append([e["id"], s["query"], v["gene"], v["hgvs_c_full"] or v["c"], v["p"], v["type"], v["zyg"],
                             v["cls"], v["pts"], v["crit"], s["status"]])
    with open(os.path.join(folder, "cohort_summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["patient_id", "sequence", "gene", "hgvs_c", "hgvs_p", "type", "zygosity", "classification",
                    "points", "criteria", "patient_result"])
        w.writerows(rows)
    trs = []
    for e, ps in prow:
        if ps is None:
            trs.append(f"<tr><td><b>{escape(e['id'])}</b></td><td colspan=7 style='color:#b45309'>Analysis failed: "
                       f"{escape(e.get('error', ''))}</td></tr>")
            continue
        pt = ps["patient"]
        vs = list({(v["c"], v["p"]): v for s in ps["samples"] for v in s["variants"]}.values())
        vtxt = "<br>".join(
            f"<span class='mono'>{escape(v['c'])} {escape(v['p'])}</span> "
            f"<span class='badge' style='background:{CLASS_COLORS.get(v['cls'], '#999')}'>{escape(v['cls'])}</span>"
            for v in vs) or "No variants detected"
        status = "; ".join(sorted({s["status"] for s in ps["samples"]}))
        fam_txt = (f"{escape(e.get('family', ''))}<br><span class='muted'>{escape(e.get('relationship', ''))}"
                   f"{', ' + escape(e.get('affected', '')) if e.get('affected') else ''}</span>") if e.get("family") else ""
        trs.append(f"<tr><td><b><a href='{escape(e['folder'])}/report.html'>{escape(e['id'])}</a></b></td><td>{fam_txt}</td>"
                   f"<td>{escape(pt.get('sex', ''))}</td><td>{escape(pt.get('age', ''))}</td>"
                   f"<td>{escape(pt.get('indication', ''))}<br><span class='muted'>{escape('; '.join(ps['phenotypes']))}</span></td>"
                   f"<td>{vtxt}</td><td>{escape(status)}</td></tr>")
    fam_html = ""
    if summary.get("families"):
        from .family import family_html
        for fa in summary["families"]:
            fa2 = dict(fa, evidence={})
            for e in fa["evidence"]:
                fa2["evidence"].setdefault(e["key"], []).append((e["code"], e["strength"], e["rationale"]))
            fam_html += family_html(fa2)
        fam_html = "<h2>Family and couple analysis</h2>" + fam_html
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>Cohort summary - {escape(summary['name'])}</title>
<style>body{{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:24px;color:#1d2433}}
table{{border-collapse:collapse;width:100%;font-size:13px}}td,th{{border-bottom:1px solid #e5e7eb;padding:7px;text-align:left;vertical-align:top}}
th{{background:#f3f4f6}}.badge{{color:#fff;padding:1px 7px;border-radius:9px;font-size:11px}}.mono{{font-family:Menlo,Consolas,monospace}}
.muted{{color:#6b7280;font-size:12px}}
.family{{border:1px solid #e5e7eb;border-radius:10px;padding:10px 16px;margin:16px 0}}li.risk{{color:#b2182b}}
.warn{{color:#b45309}}</style></head><body>
<h2>Cohort summary</h2><p class="muted">CLAVAGEN v{escape(summary['version'])} | {escape(summary['date'])} |
{len(summary['patients'])} patient(s) | online annotation: {'on' if summary['online'] else 'off'}</p>
<table><tr><th>Patient</th><th>Family</th><th>Sex</th><th>Age</th><th>Indication / phenotype</th><th>Variants (ACMG class)</th><th>Result</th></tr>
{''.join(trs)}</table>
{fam_html}<p class="muted">Click a patient ID to open the full report. Decision support only - all
classifications require expert review.</p></body></html>"""
    with open(os.path.join(folder, "cohort_summary.html"), "w") as fh:
        fh.write(html)
