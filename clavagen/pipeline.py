"""End-to-end pipeline: inputs -> alignment -> variants -> annotation -> ACMG -> outputs."""

from __future__ import annotations

import datetime as _dt
import platform

from . import __version__
from .acmg import CLASS_ABBR, combine, evaluate, parse_user_evidence
from .align import ScoringScheme, align
from .annotation import context_for_reference, load_annotation
from .caller import call_variants
from .consequence import CONSEQUENCE_LABEL, CONSEQUENCE_SEVERITY, TranscriptModel, annotate_dna, annotate_protein
from .genetic_code import AA_CLASS, HYDROPATHY, blosum62, grantham, grantham_category
from .knowledge import ACMG_SF_VERSION
from .seqio import load_sequences

CLASS_ORDER = ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]


def _parse_range(s):
    if not s:
        return None
    if isinstance(s, (list, tuple)):
        return int(s[0]) - 1, int(s[1])
    a, b = s.replace("..", "-").split("-")
    return int(a) - 1, int(b)


def build_transcript(ref_seq, ctx, cds_override=None):
    cds = _parse_range(cds_override) or _parse_range(ctx.get("cds"))
    exons = [(int(s) - 1, int(e)) for s, e in ctx.get("exons", [])] or None
    return TranscriptModel(ref_seq, cds=cds, exons=exons, exon_boundaries_c=ctx.get("exon_boundaries_c"))


def _ortholog_columns(ref_prot, orthologs):
    """Map each reference residue index -> list of ortholog residues (global alignment)."""
    cols = {i: [] for i in range(len(ref_prot))}
    for o in orthologs:
        aln = align("ref", ref_prot, o.id, o.seq.rstrip("*"), mol_type="protein", mode="global")
        ri = 0
        for a, b in zip(aln.aligned_ref, aln.aligned_query):
            if a != "-":
                if b != "-":
                    cols[ri].append(b)
                ri += 1
    return cols


class Pipeline:
    def __init__(self, queries, references, annotation=None, gene=None, cds=None, pairs=None,
                 orthologs=None, evidence=None, phenotypes=None, mol_type="auto", aligner="auto",
                 mode="local", strict=False, patient=None, min_identity=80.0, online=False,
                 known_variants=None, log=None, auto_cache=None, use_clinvar=True, clinvar_min_stars=0,
                 family_evidence=None, exclude_known=None, default_reference=None, phases=None,
                 clinical_evidence=None):
        ft = None if mol_type == "auto" else mol_type
        self.log = log or (lambda *_: None)
        self.online = online
        self.queries = load_sequences(queries, "query", ft)
        self.references = self._load_references(references)
        if not self.queries:
            raise ValueError("No query sequences supplied.")
        if not self.references:
            raise ValueError("No reference sequences supplied.")
        self.annotation_path = annotation
        self.ann = load_annotation(annotation) if isinstance(annotation, str) else (annotation or {})
        self.gene = gene
        self.cds = cds
        self.pairs = dict(p.split(":", 1) for p in (pairs or []))
        self.orthologs = load_sequences(orthologs, "ortholog", "protein") if orthologs else []
        self.user_ev = parse_user_evidence(evidence)
        # guided clinical evidence: {variant key or "*": {...}} (see clavagen.evidence); converted per gene
        # rule profile when the reference's annotation has one
        self.clinical_raw = clinical_evidence or {}
        self.phenotypes = [p.strip() for p in (phenotypes or []) if p.strip()]
        self.aligner = aligner
        self.mode = mode
        self.strict = strict
        self.patient = patient or {}
        self.min_identity = min_identity
        self.manual_known = list(known_variants or [])
        self._auto_cache = auto_cache if auto_cache is not None else {}
        self.use_clinvar = use_clinvar
        self.family_evidence = family_evidence or {}
        self.phases = phases or []
        self.exclude_known = set(exclude_known or [])
        self.default_reference = (default_reference or "").strip() or None
        self.clinvar_min_stars = clinvar_min_stars

    def _load_references(self, specs):
        """References may be files/literals, or RefSeq accessions (e.g. NM_000518.5) fetched online."""
        import os
        import re as _re
        from .seqio import SeqRecord
        import hashlib
        out, rest = [], []
        for sp in specs or []:
            s = sp.strip()
            if not os.path.isfile(s) and _re.fullmatch(r"N[MR]_\d+(\.\d+)?", s):
                if not self.online:
                    raise ValueError(f"Reference {s} is an accession: enable online annotation to download it.")
                from .online import fetch_refseq
                self.log(f"Downloading reference {s} from NCBI RefSeq...")
                gb = fetch_refseq(s)
                rec = SeqRecord(id=gb.get("transcript") or s, seq=gb["seq"],
                                description=f"({gb.get('gene')}) RefSeq {gb.get('transcript')} downloaded from NCBI",
                                mol_type="dna", source=f"NCBI nuccore {s}",
                                sha256=hashlib.sha256(gb["seq"].encode()).hexdigest())
                rec.genbank = gb
                out.append(rec)
            else:
                rest.append(sp)
        return out + load_sequences(rest, "ref", None)

    def _context(self, r):
        """Annotation context for reference r, including online auto-annotation."""
        auto = None
        online_info = None
        if self.online:
            from .annotation import user_context
            from .online import auto_annotation, detect_gene_and_transcript
            uc = user_context(self.ann, r.id)
            if self.gene:
                uc = dict(uc, gene=self.gene)
            gb = getattr(r, "genbank", None)
            if gb:
                uc.setdefault("gene", gb.get("gene"))
                uc.setdefault("transcript", gb.get("transcript"))
            gene, tx = detect_gene_and_transcript(r, uc, resolve=True, log=self.log)
            if gene:
                key = (gene, tx)
                if key not in self._auto_cache:
                    self.log(f"Searching online databases for {gene}{' / ' + tx if tx else ''}...")
                    self._auto_cache[key] = auto_annotation(gene, tx, log=self.log)
                auto = self._auto_cache[key]
                online_info = auto.get("_online")
            else:
                online_info = {"sources": [], "errors": ["Gene could not be identified - give a gene symbol, "
                                                         "an NM_ accession in the reference header, or an annotation file."]}
        ctx = context_for_reference(self.ann, r.id, self.gene, auto=auto, manual_known=self.manual_known)
        if not ctx.get("gene"):
            # offline (or nothing found online): still take the gene symbol from the reference header
            # ('gene=CELF2', '(CELF2)', 'GN=CELF2' or a built-in gene name), as the window does
            from .annotation import user_context
            from .online import detect_gene_and_transcript
            hg, _ = detect_gene_and_transcript(r, user_context(self.ann, r.id), resolve=False)
            if hg:
                ctx = context_for_reference(self.ann, r.id, hg, auto=auto, manual_known=self.manual_known)
        if self.exclude_known:  # leave-one-out (benchmarking)
            from .annotation import gene_info, normalize_hgvs, user_context
            from .online import recalc_gene_evidence
            # drop the tested variant's own record(s): by DNA change, and records that give only the protein
            # change (they could otherwise match the variant as "same amino-acid change", i.e. PS1 for itself)
            ctx["known_variants"] = [
                k for k in ctx.get("known_variants", [])
                if normalize_hgvs(k.get("hgvs_c", "")) not in self.exclude_known
                and not (not k.get("hgvs_c") and normalize_hgvs(k.get("hgvs_p", "")) in self.exclude_known)]
            # gene-level evidence derived from ClinVar (hotspots, LOF mechanism) is recounted without the
            # tested variant; values from the built-in gene list or the user's annotation are kept
            fixed_hs, fixed_lof = set(), None
            for src in (gene_info(ctx.get("gene")) if ctx.get("gene") else {}, user_context(self.ann, r.id)):
                fixed_hs |= set((src or {}).get("hotspots") or [])
                if (src or {}).get("lof_mechanism") is not None:
                    fixed_lof = src["lof_mechanism"]
            recalc_gene_evidence(ctx, fixed_hs, fixed_lof)
        # CDS/exon structure from RefSeq if the user did not define it
        gb = getattr(r, "genbank", None) or (auto or {}).get("_refseq")
        if gb and gb.get("cds") and not ctx.get("cds") and not self.cds and r.mol_type == "dna":
            cs, ce = gb["cds"]
            cds_seq = gb["seq"][cs - 1:ce]
            if r.seq == gb.get("seq"):
                ctx["cds"] = [cs, ce]
            else:
                i = r.seq.find(cds_seq)
                if i >= 0:
                    ctx["cds"] = [i + 1, i + len(cds_seq)]
            if ctx.get("cds") and not ctx.get("exons") and gb.get("exon_boundaries_c") is not None:
                ctx.setdefault("exon_boundaries_c", gb["exon_boundaries_c"])
        if auto and auto.get("_critical_sites"):
            ctx["critical_sites"] = auto["_critical_sites"]
        import copy
        ctx["_online"] = copy.deepcopy(online_info)
        ctx["_derived"] = (auto or {}).get("_derived", [])
        return ctx

    def _resolve_het_indel(self, trace, q, r, ref_seq, aln, scheme, notes):
        """Split a trace with a heterozygous indel into two alleles and call the second allele."""
        from .caller import call_variants
        # map trimmed trace index -> reference base using the current alignment
        qmap = {}
        rp, qp = aln.ref_start, aln.query_start
        for a, b in zip(aln.aligned_ref, aln.aligned_query):
            if a != "-" and b != "-":
                qmap[qp] = ref_seq[rp]
            if a != "-":
                rp += 1
            if b != "-":
                qp += 1
        from .trace import deconvolute
        a1, a2 = deconvolute(trace, lambda i: qmap.get(i))
        q.seq = a1
        aln1 = align(r.id, ref_seq, q.id, a1, mol_type="dna", mode=self.mode, scheme=scheme, backend=self.aligner)
        aln2 = align(r.id, ref_seq, q.id + "_allele2", a2, mol_type="dna", mode=self.mode, scheme=scheme,
                     backend=self.aligner)
        raw2, _ = call_variants(aln2, ref_seq, "dna")
        # reference position where the double peaks start
        start_ref = None
        rp, qp = aln1.ref_start, aln1.query_start
        for a, b in zip(aln1.aligned_ref, aln1.aligned_query):
            if b != "-" and qp == trace.het_indel_site:
                start_ref = rp
                break
            if a != "-":
                rp += 1
            if b != "-":
                qp += 1
        extra = []
        for v in raw2:
            if v.kind in ("deletion", "insertion", "duplication", "delins") and \
                    (start_ref is None or abs(v.ref_start - start_ref) <= 40):
                v.zygosity = "heterozygous"
                v.notes.append("Heterozygous indel resolved from double peaks in the Sanger trace (second allele "
                               "reconstructed from secondary peaks).")
                extra.append(v)
        if extra:
            notes.append(f"Sanger trace: heterozygous indel detected from sustained double peaks starting at trace "
                         f"position {trace.trim[0] + trace.het_indel_site + 1}; alleles separated "
                         f"({', '.join(f'{x.kind} at ref {x.ref_start + 1}' for x in extra)}).")
        else:
            notes.append("Sanger trace: sustained double peaks detected but no indel could be resolved - check for "
                         "a mixed sample or poor-quality read.")
        trace.allele2 = a2
        return aln1, extra

    def _trace_views(self, trace, aln, ref_seq, variants):
        """Attach a chromatogram SVG around each variant (DNA input from a Sanger trace)."""
        from .trace import chromatogram_svg
        r2q = {}
        rp, qp = aln.ref_start, aln.query_start
        for a, b in zip(aln.aligned_ref, aln.aligned_query):
            if a != "-":
                r2q[rp] = qp
                rp += 1
            if b != "-":
                qp += 1
        for v in variants:
            r0 = v["ref_pos"] - 1
            qi = r2q.get(r0)
            if qi is None:
                continue
            ref_win = "".join(ref_seq[r0 + j] if 0 <= r0 + j < len(ref_seq) else " " for j in range(-12, 13))
            v["trace_svg"] = chromatogram_svg(trace, qi, ref_win,
                                              title=f"{v.get('hgvs_c_short') or v.get('hgvs_ref')} - "
                                                    f"{v['zygosity'].split(' ')[0]}")

    def _infer_nucleotide(self, ctx, variants, ref_prot, notes, tm=None):
        """Protein input: derive the nucleotide change(s) that explain each amino-acid change from the
        transcript (DNA reference or RefSeq), so population / in-silico / ClinVar lookups can use them.
        Single-nucleotide changes for missense/nonsense; small insertions/deletions for frameshifts."""
        from .genetic_code import CODON_TABLE, translate
        from .truncation import infer_frameshift_dna
        gb = ctx.get("_refseq")
        if tm is not None and getattr(tm, "cds_seq", ""):
            cds, utr3 = tm.cds_seq, getattr(tm, "utr3_seq", "")
            src = ctx.get("transcript") or "the reference transcript"
        elif gb and gb.get("cds") and gb.get("seq"):
            cs, ce = gb["cds"]
            cds, utr3 = gb["seq"][cs - 1:ce], gb["seq"][ce:]
            src = gb.get("transcript")
        else:
            notes.append("Protein input: no RefSeq transcript available, so the nucleotide change cannot be inferred "
                         "(population frequency and REVEL/SpliceAI lookups need it).")
            return
        tprot = translate(cds).rstrip("*")
        off = tprot.find(ref_prot.rstrip("*"))
        if off < 0:
            notes.append(f"Protein input: reference protein does not match the transcript {src} "
                         "translation; nucleotide inference skipped.")
            return
        tx = ctx.get("transcript") or (gb or {}).get("transcript")
        n_inf = 0
        for v in variants:
            if v.get("kind") == "frameshift" and v.get("aa_pos"):
                cands = infer_frameshift_dna(cds, utr3, off, v)
                if len(cands) == 1:
                    v["hgvs_c"] = v["hgvs_c_short"] = cands[0]
                    v["hgvs_c_full"] = f"{tx}:{cands[0]}" if tx else cands[0]
                    v["nucleotide_inferred"] = True
                    v["notes"].append(f"Frameshift explained by {cands[0]}: the only small insertion/deletion "
                                      f"(1-5 nt) near codon {v['aa_pos']} that reproduces the observed new "
                                      "residues; confirm by DNA sequencing.")
                    n_inf += 1
                elif cands:
                    v["hgvs_c_candidates"] = cands
                    v["notes"].append(f"Frameshift could be explained by: {', '.join(cands)} (DNA needed to decide).")
                else:
                    v["notes"].append("No 1-5 nt insertion/deletion near the break reproduces the new residues; "
                                      "the change may be complex (e.g. splicing or a larger rearrangement).")
                continue
            if v.get("kind") != "SNV" or not v.get("aa_pos") or v.get("consequence") not in ("missense_variant", "stop_gained"):
                continue
            pos = v["aa_pos"] + off
            codon = cds[(pos - 1) * 3:pos * 3]
            if len(codon) < 3:
                continue
            alt = v.get("alt_aa") or ("*" if v["consequence"] == "stop_gained" else "")
            cands = []
            for k in range(3):
                for b in "ACGT":
                    if b == codon[k]:
                        continue
                    newc = codon[:k] + b + codon[k + 1:]
                    if CODON_TABLE.get(newc) == alt:
                        cands.append(f"c.{(pos - 1) * 3 + k + 1}{codon[k]}>{b}")
            v["ref_codon"] = codon
            if len(cands) == 1:
                v["hgvs_c"] = v["hgvs_c_short"] = cands[0]
                v["hgvs_c_full"] = f"{tx}:{cands[0]}" if tx else cands[0]
                v["nucleotide_inferred"] = True
                newc = list(codon)
                k = int(cands[0][2:-3]) - (pos - 1) * 3 - 1
                newc[k] = cands[0][-1]
                v["codon_change"] = f"{codon}>{''.join(newc)}"
                v["notes"].append(f"Nucleotide change {cands[0]} inferred from protein (only single-nucleotide change "
                                  f"from reference codon {codon} giving {v['hgvs_p']}); confirm by DNA sequencing.")
                n_inf += 1
            elif cands:
                v["hgvs_c_candidates"] = cands
                v["notes"].append(f"Possible nucleotide changes for {v['hgvs_p']}: {', '.join(cands)} (DNA needed to decide).")
            else:
                v["notes"].append(f"{v['hgvs_p']} needs more than one nucleotide change from codon {codon}.")
        if n_inf:
            notes.append(f"Protein input: nucleotide change inferred for {n_inf} variant(s) from {tx} "
                         "(single-nucleotide explanation). Zygosity cannot be determined from protein.")

    def _vep(self, ctx, variants, notes):
        """Per-variant gnomAD AF / REVEL / SpliceAI / AlphaMissense via Ensembl VEP."""
        from .annotation import normalize_hgvs
        from .online import mappable, transcript_mapping, vep_lookup
        tx = ctx.get("transcript")
        if not tx:
            notes.append("Variant-level online lookup skipped: transcript accession unknown.")
            return
        tmap = transcript_mapping(tx)
        if tmap["note"]:
            notes.append(tmap["note"])

        def acc(c):  # current version when the c. position carries over unchanged (see transcript_mapping)
            return tmap["acc"] if tmap["acc"] != tx and mappable(c, tmap) else tx
        todo = {f"{acc(v['hgvs_c'])}:{v['hgvs_c']}": v for v in variants if v.get("hgvs_c")}
        multi = {}
        for v in variants:
            if not v.get("hgvs_c") and v.get("hgvs_c_candidates"):
                for c in v["hgvs_c_candidates"]:
                    multi[f"{acc(c)}:{c}"] = v
        if not todo and not multi:
            return
        todo_all = dict(todo, **multi)
        unresolved = ctx.setdefault("population_unresolved", set())
        try:
            self.log(f"Looking up {len(todo)} variant(s) in gnomAD/REVEL/SpliceAI via Ensembl VEP...")
            res = vep_lookup(list(todo_all), tx)
        except Exception as e:
            res = {}
            err = e
        else:
            err = None
        # VEP's RefSeq set can lag behind NCBI (e.g. NM_000249.5 at NCBI, NM_000249.4 in VEP): retry what failed
        # on the other version before giving up
        retry = {}
        for h, v in todo_all.items():
            if h not in res and tmap["acc"] != tx:
                acc_used, c = h.split(":", 1)
                other = tx if acc_used == tmap["acc"] else (tmap["acc"] if mappable(c, tmap) else None)
                if other:
                    retry[f"{other}:{c}"] = h
        if retry:
            try:
                res2 = vep_lookup(list(retry), tx)
            except Exception as e:
                res2, err = {}, err or e
            for h2, h in retry.items():
                if h2 in res2:
                    res[h] = res2[h2]
                    res[h]["_acc"] = h2.split(":", 1)[0]
        for h, v in todo.items():
            v["vep_resolved"] = h in res
            v["lookup_transcript"] = (res[h].get("_acc") or h.split(":", 1)[0]) if h in res else ""
            if h not in res:
                unresolved.add(normalize_hgvs(v["hgvs_c"]))
        if err is not None and not res:
            notes.append(f"Ensembl VEP lookup failed ({err}); population/in-silico data not added.")
            if ctx.get("_online"):
                ctx["_online"].setdefault("errors", []).append(f"Ensembl VEP: {err}")
            return
        ctx.setdefault("population", {})
        ctx.setdefault("insilico", {})
        n_af = 0
        for h, v in todo.items():  # noqa
            r = res.get(h)
            if r is None:
                continue
            key = normalize_hgvs(v["hgvs_c"])
            if "af" in r and key not in ctx["population"]:
                ctx["population"][key] = {"af": r["af"], "source": "gnomAD (Ensembl VEP)", "detail": r.get("af_detail")}
                n_af += 1
            ins = {k: r[k] for k in ("REVEL", "SpliceAI", "AlphaMissense") if k in r}
            if ins and key not in ctx["insilico"]:
                ctx["insilico"][key] = ins
            if r.get("clinvar_vep"):
                v.setdefault("database_assertions", []).append(f"ClinVar (via VEP): {r['clinvar_vep']}")
            if r.get("rsid"):
                v["rsid"] = r["rsid"]
        # protein input with several possible nucleotide changes: aggregate under the p. key
        agg = {}
        for h, v in multi.items():
            r = res.get(h) or {}
            a = agg.setdefault(id(v), {"v": v, "af": None, "ins": {}})
            if "af" in r:
                a["af"] = max(a["af"] or 0, r["af"])
            for k in ("REVEL", "SpliceAI", "AlphaMissense"):
                if k in r:
                    a["ins"][k] = max(a["ins"].get(k, 0), r[k])
        for a in agg.values():
            key = normalize_hgvs(a["v"]["hgvs_p"])
            if a["af"] is not None:
                ctx["population"].setdefault(key, {"af": a["af"], "source": "gnomAD (Ensembl VEP; max over possible "
                                                                           "nucleotide changes)"})
            if a["ins"]:
                ctx["insilico"].setdefault(key, a["ins"])
        if not res:
            notes.append("Ensembl VEP could not resolve the variant(s) on this transcript; "
                         "population/in-silico data not added.")
        elif not ctx.get("population_dataset"):
            ctx["population_dataset"] = "gnomAD v4 exomes/genomes (via Ensembl VEP)"
        if ctx.get("_online") is not None:
            ctx["_online"].setdefault("sources", []).append(
                f"Ensembl VEP: {len(res)}/{len(todo)} variants annotated, {n_af} with gnomAD AF")

    # ------------------------------------------------------------ pairing
    def _find_ref(self, rid):
        base = rid.split(".")[0]
        for r in self.references:
            if r.id == rid or r.id.split(".")[0] == base:
                return r
        for r in self.references:  # gene symbol given instead of an ID
            if f"({rid.upper()})" in f"{r.id} {r.description}".upper():
                return r
        return None

    def _choose_reference(self, q):
        if self.default_reference and q.id not in self.pairs:
            r = self._find_ref(self.default_reference)
            if r is None:
                raise ValueError(f"Reference '{self.default_reference}' chosen for this patient is not loaded.")
            return r, None
        if q.id in self.pairs:
            for r in self.references:
                if r.id == self.pairs[q.id]:
                    return r, None
            raise ValueError(f"Reference '{self.pairs[q.id]}' for query '{q.id}' not found.")
        if len(self.references) == 1:
            return self.references[0], None
        best = None
        for r in self.references:
            rseq, mol = self._ref_sequence_for(q, r)
            if rseq is None:
                continue
            aln = align(r.id, rseq, q.id, q.seq, mol_type=mol, mode="local", backend=self.aligner)
            if best is None or aln.score > best[1].score:
                best = (r, aln)
        if best is None:
            raise ValueError(f"No compatible reference for query '{q.id}'.")
        return best

    def _ref_sequence_for(self, q, r):
        if q.mol_type == r.mol_type:
            return r.seq, q.mol_type
        if q.mol_type == "protein" and r.mol_type == "dna":
            tm = build_transcript(r.seq, self._context(r), self.cds)
            return tm.protein.rstrip("*"), "protein"
        return None, None

    # ------------------------------------------------------------ run
    def run(self):
        samples = []
        for qi, q in enumerate(self.queries, 1):
            r, pre = self._choose_reference(q)
            ctx = self._context(r)
            tm = None
            notes = list(q.warnings)
            if r.mol_type == "dna":
                tm = build_transcript(r.seq, ctx, self.cds)
                notes += tm.warnings
            if q.mol_type == "protein":
                mode = "protein"
                ref_seq = tm.protein.rstrip("*") if (r.mol_type == "dna" and tm) else r.seq.rstrip("*")
                if r.mol_type == "dna":
                    notes.append("Protein query compared with the translated reference CDS.")
            else:
                mode = "dna"
                ref_seq = r.seq
                if r.mol_type == "protein":
                    raise ValueError(f"DNA query '{q.id}' cannot be compared with protein reference '{r.id}'.")
            scheme = ScoringScheme.default("protein" if mode == "protein" else "dna")
            aln = pre if (pre is not None and pre.aligned_ref and self.mode == "local" and
                          pre.mol_type == scheme.mol_type) else \
                align(r.id, ref_seq, q.id, q.seq, mol_type=scheme.mol_type, mode=self.mode, scheme=scheme,
                      backend=self.aligner)
            trace = getattr(q, "trace", None)
            if mode == "dna":
                # reverse reads: try the reverse complement and keep the better alignment
                from .genetic_code import reverse_complement
                aln_rc = align(r.id, ref_seq, q.id, reverse_complement(q.seq), mol_type="dna", mode=self.mode,
                               scheme=scheme, backend=self.aligner)
                if aln_rc.score > aln.score * 1.2:
                    aln = aln_rc
                    if trace is not None:
                        trace.reverse_complement()
                        q.seq = trace.consensus()
                        aln = align(r.id, ref_seq, q.id, q.seq, mol_type="dna", mode=self.mode, scheme=scheme,
                                    backend=self.aligner)
                    else:
                        q.seq = reverse_complement(q.seq)
                    notes.append("Query is in reverse orientation relative to the reference - reverse complement used.")
            extra = []
            if trace is not None and trace.het_indel_start is not None and mode == "dna":
                aln, extra = self._resolve_het_indel(trace, q, r, ref_seq, aln, scheme, notes)
            raw, qc = call_variants(aln, ref_seq, scheme.mol_type, merge_adjacent=(mode == "dna"))
            trunc_warn = []
            if mode == "protein":
                from .truncation import detect_protein_truncation
                tv, tnotes, trunc_warn = detect_protein_truncation(
                    aln, ref_seq, q.seq,
                    align_fn=lambda rs, qs: align("ref", rs, "tail", qs, mol_type="protein", mode="local",
                                                  backend=self.aligner))
                if tv is not None and tv.kind == "frameshift":
                    raw = [x for x in raw if x.ref_start < tv.ref_start]   # absorbed into the frameshift
                if tv is not None and not any(x.ref_start == tv.ref_start for x in raw):
                    tv.notes.extend(tnotes)
                    raw.append(tv)
                    raw.sort(key=lambda x: x.ref_start)
            if extra:
                have = {(x.ref_start, x.ref, x.alt) for x in raw}
                raw += [x for x in extra if (x.ref_start, x.ref, x.alt) not in have]
                raw.sort(key=lambda x: x.ref_start)
            ortho = None
            ref_prot = (tm.protein.rstrip("*") if tm else "") if mode == "dna" else ref_seq
            if self.orthologs and ref_prot:
                ortho = _ortholog_columns(ref_prot, self.orthologs)
            variants = []
            for vi, rv in enumerate(raw, 1):
                v = self._molecular(rv, mode, tm, ref_seq, ctx, ortho)
                v["id"] = f"Q{qi}V{vi}"
                v["query_id"], v["ref_id"] = q.id, r.id
                variants.append(v)
            if mode == "protein" and variants:
                self._infer_nucleotide(ctx, variants, ref_seq, notes, tm)
            if self.online and variants:
                self._vep(ctx, variants, notes)
            for v in variants:
                self._classify(v, ctx)
            on = ctx.get("_online")
            if on:
                for e in on.get("errors", []):
                    notes.append(f"Online annotation: {e}")
            self._qc(q, aln, tm, variants, mode, qc)
            qc["warnings"].extend(trunc_warn)
            qc["notes"] = notes + qc.get("notes", [])
            if trace is not None:
                self._trace_views(trace, aln, ref_seq, variants)
            sample = {"query": q, "reference": r, "ref_seq": ref_seq, "mode": mode, "alignment": aln,
                      "scheme": scheme, "tm": tm, "context": ctx, "variants": variants, "qc": qc}
            sample["clinical"] = clinical_interpretation(sample, self.patient)
            samples.append(sample)
        from .genotype import apply_patient_genotype
        apply_patient_genotype(samples, self.patient, self.phases)
        return {
            "meta": {
                "software": "CLAVAGEN", "version": __version__,
                "run_date": _dt.datetime.now().isoformat(timespec="seconds"),
                "python": platform.python_version(), "aligner_mode": self.mode,
                "strict_mode": self.strict, "annotation_file": self.annotation_path,
                "online_annotation": self.online, "manual_known_variants": len(self.manual_known),
                "acmg_sf_version": ACMG_SF_VERSION, "patient": self.patient,
                "phenotypes": self.phenotypes,
                "inputs": [{"id": s.id, "role": "query", "type": s.mol_type, "length": len(s), "sha256": s.sha256,
                            "source": s.source} for s in self.queries] +
                          [{"id": s.id, "role": "reference", "type": s.mol_type, "length": len(s),
                            "sha256": s.sha256, "source": s.source} for s in self.references],
                "guidelines": ["ACMG/AMP 2015 (Richards et al.)", "ClinGen SVI: PVS1 tree, PM2_Supporting, "
                               "PP5/BP6 retired, REVEL calibration (Pejaver 2022)", "Tavtigian points (2020)"],
            },
            "samples": samples,
        }

    # ------------------------------------------------------------ helpers
    def _molecular(self, rv, mode, tm, ref_seq, ctx, ortho):
        if mode == "dna":
            a = annotate_dna(rv, tm)
        else:
            a = annotate_protein(rv, ref_seq)
        v = {
            "mol_type": mode, "kind": rv.kind, "ref_pos": rv.ref_start + 1,
            "ref_end": rv.ref_end if rv.ref_end > rv.ref_start else rv.ref_start,
            "ref_allele": rv.ref, "alt_allele": rv.alt, "zygosity": rv.zygosity,
            "repeat_unit": rv.repeat_unit, "repeat_copies": rv.repeat_copies, "notes": list(rv.notes),
            "query_pos": rv.query_pos + 1,
        }
        v.update(a)
        if mode == "protein":
            v["hgvs_c"] = ""
            v["zygosity"] = "not determinable (protein input)"
        tx = ctx.get("transcript")
        pid = ctx.get("protein_id")
        v["hgvs_c_short"] = v.get("hgvs_c", "")
        v["hgvs_c_full"] = f"{tx}:{v['hgvs_c']}" if tx and v.get("hgvs_c") else v.get("hgvs_c", "")
        v["hgvs_p_full"] = f"{pid}:{v['hgvs_p']}" if pid and v.get("hgvs_p") else v.get("hgvs_p", "")
        v["gene"] = ctx.get("gene", "")
        z = v["zygosity"]
        v["zygosity_short"] = "het" if z == "heterozygous" else ("hom (apparent)" if z.startswith("homozygous") else "n.d.")
        v["consequence_label"] = CONSEQUENCE_LABEL.get(v["consequence"], v["consequence"])
        v["junctions_c"] = tm.junctions_c if tm else list(ctx.get("exon_boundaries_c") or [])
        # amino-acid level metrics
        if v["consequence"] == "missense_variant":
            r, al = v["ref_aa"], v["alt_aa"]
            v["blosum62"] = blosum62(r, al)
            g = grantham(r, al)
            v["grantham"] = g
            v["grantham_category"] = grantham_category(g)
            v["aa_property_change"] = f"{AA_CLASS.get(r, '?')} -> {AA_CLASS.get(al, '?')}"
            if r in HYDROPATHY and al in HYDROPATHY:
                v["hydropathy_change"] = round(HYDROPATHY[al] - HYDROPATHY[r], 2)
        if ortho and v.get("aa_pos") and (v["aa_pos"] - 1) in ortho:
            col = ortho[v["aa_pos"] - 1]
            ra = v.get("ref_aa", "")
            v["conservation"] = {
                "n_orthologs": len(col),
                "ref_identity": round(sum(1 for c in col if c == ra) / len(col), 3) if col else 0,
                "alt_observed": bool(v.get("alt_aa")) and v.get("alt_aa") in col,
                "column": "".join(col),
            }
        pos = v.get("aa_pos")
        doms = ctx.get("domains") or []
        v["domains"] = [d for d in doms if pos and int(d["start"]) <= pos <= int(d["end"])]
        ptc = v.get("ptc_pos")
        v["domains_downstream"] = [d for d in doms if ptc and int(d["end"]) >= ptc]
        return v

    def _classify(self, v, ctx):
        from .annotation import population_lookup, insilico_lookup
        pop = population_lookup(ctx, [v.get("hgvs_c"), v.get("hgvs_p")])
        v["population_af"] = pop.get("af") if pop else None
        v["population_source"] = (pop.get("source") if pop else None) or ctx.get("population_dataset")
        v["insilico"] = insilico_lookup(ctx, [v.get("hgvs_c"), v.get("hgvs_p")]) or {}
        from .family import vkey
        from .evidence import from_json
        from .annotation import normalize_hgvs as _nh
        ce = {(_nh(k) if k != "*" else "*"): x for k, x in from_json(self.clinical_raw, ctx.get("rules")).items()}
        crit = evaluate(v, ctx, self.user_ev, strict=self.strict, patient_phenotypes=self.phenotypes,
                        family_ev=self.family_evidence.get(vkey(v)), clinical_ev=ce)
        cls = combine(crit)
        v["criteria"] = [c.to_dict() for c in crit]
        v.update(cls)
        v["acmg_classification"] = cls["classification"]
        v["classification_basis"] = "ACMG/AMP criteria (automated)"
        cm = v.get("clinvar_match")
        v["clinvar_classification"] = cm.get("classification") if cm else ""
        v["clinvar_stars"] = cm.get("stars") if cm else None
        v["clinvar_review"] = cm.get("review_status", "") if cm else ""
        v["clinvar_source"] = cm.get("source", "") if cm else ""
        v["clinvar_conditions"] = cm.get("conditions", []) if cm else []
        dbc = normalize_db_class(v["clinvar_classification"])
        is_clinvar = str(cm.get("source", "")).startswith("ClinVar") if cm else False
        stars = (cm.get("stars") or 0) if cm else 0
        if cm and not dbc and v["clinvar_classification"]:
            v["notes"].append(f"Database classification '{v['clinvar_classification']}' is conflicting/unclear - "
                              "automated ACMG class reported; review the submissions.")
        if (self.use_clinvar and cm and dbc and (stars >= self.clinvar_min_stars or not is_clinvar)
                and "conflicting" not in (cm.get("review_status") or "").lower()):
            v["classification"] = dbc
            if is_clinvar:
                level = {0: "no assertion criteria - weak", 1: "limited review"}.get(stars, "reviewed")
                v["classification_basis"] = (f"{cm.get('source', 'ClinVar')} {stars}\u2605 ({cm.get('review_status') or 'review status n/a'}; "
                                             f"{level})")
            else:
                v["classification_basis"] = f"Known-variant list: {cm.get('source', 'user list')}"
            if cm.get("match") != "exact":
                v["classification_basis"] += " - matched on the protein change; nucleotide not confirmed"
            v["clinvar_review_level"] = level if is_clinvar else "user list"
            if dbc != cls["classification"]:
                v["notes"].append(f"Reported class taken from reviewed ClinVar record ({dbc}); automated ACMG criteria "
                                  f"alone give {cls['classification']} ({cls['points']} points) because literature "
                                  "evidence (functional, case and segregation data) is not available to the software.")
            v["vus_subtier"] = "" if dbc != "Uncertain significance" else v.get("vus_subtier", "")
        v["class_abbr"] = CLASS_ABBR.get(v["classification"], "?")
        v.setdefault("database_assertions", [])
        v["interpretation"] = variant_sentence(v, ctx)
        return v

    def _qc(self, q, aln, tm, variants, mode, qc):
        qc["alignment"] = aln.stats
        qc.setdefault("warnings", [])
        st = aln.stats
        if st["identity_pct"] < self.min_identity:
            qc["warnings"].append(f"Low identity ({st['identity_pct']}%) - wrong reference, wrong strand or poor-quality sequence?")
        if mode == "dna" and tm and tm.coding:
            covered = sum(1 for g in tm.coding if aln.ref_start <= g < aln.ref_end)
            qc["cds_coverage_pct"] = round(100 * covered / len(tm.coding), 1)
            if qc["cds_coverage_pct"] < 100:
                qc["warnings"].append(f"Query covers {qc['cds_coverage_pct']}% of the CDS; variants outside the "
                                      "aligned region cannot be assessed.")
        n_n = q.seq.count("N") if mode == "dna" else q.seq.count("X")
        qc["no_call_count"] = n_n
        if aln.stats["alignment_length"] and len(variants) / aln.stats["alignment_length"] > 0.05:
            qc["warnings"].append("High variant density (>5% of aligned positions) - check sequence quality/reference.")
        if mode == "dna":
            from .seqio import IUPAC_DNA  # noqa
            qc["heterozygous_calls"] = sum(1 for v in variants if v["zygosity"] == "heterozygous")
            if getattr(q, "trace", None) is not None:
                qc.setdefault("notes", []).append(
                    "Zygosity from the Sanger trace: double peaks (secondary/primary >= 0.30) are called heterozygous; "
                    "sustained double peaks are resolved as heterozygous indels. Confirm low-level mosaicism separately.")
            else:
                qc.setdefault("notes", []).append(
                    "Zygosity is inferred from the input consensus: two-base IUPAC codes are called heterozygous; other "
                    "calls are 'apparent homozygous'. Heterozygous indels need the trace (.ab1) or NGS allele counts.")
        if qc["no_call_positions"]:
            qc["warnings"].append(f"{len(qc['no_call_positions'])} no-call (N/X) position(s) in the aligned region: "
                                  f"{', '.join(map(str, qc['no_call_positions'][:10]))}.")
        if qc.get("ambiguous_uncallable"):
            qc["warnings"].append(f"{len(qc['ambiguous_uncallable'])} ambiguous base(s) could not be resolved.")
        return qc


# ====================================================================== text
def normalize_db_class(c):
    c = (c or "").lower().strip()
    if not c or "conflicting" in c:
        return ""
    if c.startswith("pathogenic/likely pathogenic"):
        return "Likely pathogenic"
    if c.startswith("pathogenic"):
        return "Pathogenic"
    if c.startswith("likely pathogenic"):
        return "Likely pathogenic"
    if c.startswith("benign/likely benign"):
        return "Likely benign"
    if c.startswith("benign"):
        return "Benign"
    if c.startswith("likely benign"):
        return "Likely benign"
    if c.startswith("uncertain"):
        return "Uncertain significance"
    return ""


def variant_sentence(v, ctx):
    gene = ctx.get("gene") or "the reference gene"
    z = v["zygosity"]
    zyg = "Heterozygous" if z == "heterozygous" else ("Apparently homozygous" if z.startswith("homozygous") else "")
    desc = v["consequence_label"].lower()
    name = " ".join(x for x in (v.get("hgvs_c_full"), v.get("hgvs_p")) if x)
    lead = f"{zyg} {desc}" if zyg else desc.capitalize()
    s = f"{lead} variant {name} in {gene}, classified as {v['classification']}"
    if v.get("classification_basis", "").startswith("ACMG"):
        s += (f"{' (' + v['vus_subtier'] + ')' if v.get('vus_subtier') else ''} "
              f"[{', '.join(v['criteria_met']) or 'no criteria met'}; {v['points']} points].")
    else:
        s += (f" based on {v['classification_basis']}. Automated ACMG criteria: {v['acmg_classification']} "
              f"[{', '.join(v['criteria_met']) or 'no criteria met'}; {v['points']} points].")
    if v.get("clinvar_conditions"):
        s += f" ClinVar condition(s): {'; '.join(v['clinvar_conditions'][:3])}."
    if v.get("nucleotide_inferred"):
        s += " Nucleotide change inferred from protein - confirm by DNA sequencing."
    if z.startswith("not determinable"):
        s += " Zygosity not determinable from protein input."
    if v.get("database_assertions") and v.get("classification_basis", "").startswith("ACMG"):
        s += f" Previously reported: {'; '.join(v['database_assertions'])}."
    if not v.get("concordant", True):
        s += f" Note: points-based framework gives {v['points_classification']} - reviewer to adjudicate."
    return s


def clinical_interpretation(sample, patient, variants=None):
    """Interpretation for one sample. ``variants`` (V12): all of the patient's variants in this gene,
    merged across samples - so a genotype split over several files is interpreted as a whole."""
    ctx = sample["context"]
    vs = sample["variants"] if variants is None else variants
    gene = ctx.get("gene") or "(gene not specified)"
    inh = (ctx.get("inheritance") or "").upper()
    sex = (patient.get("sex") or "").lower()
    plp = [v for v in vs if v["classification"] in ("Pathogenic", "Likely pathogenic")]
    vus = [v for v in vs if v["classification"] == "Uncertain significance"]
    diseases = ctx.get("diseases") or []
    dnames = ", ".join(f"{d['name']} ({d.get('inheritance', inh or '?')})" for d in diseases) or "no disease association annotated"
    lines, recs = [], []
    status = "No reportable finding"
    if plp:
        hom = [v for v in plp if v["zygosity"].startswith("homozygous")]
        het = [v for v in plp if v["zygosity"] == "heterozygous"]
        unk = [v for v in plp if v not in hom and v not in het]
        if unk and not hom and not het:
            status = "P/LP variant identified (zygosity not determinable)"
            lines.append(f"P/LP variant(s) in {gene} detected at protein level. Zygosity cannot be determined from "
                         "protein input, and any nucleotide change shown is inferred from the reference transcript. "
                         f"Associated condition(s): {dnames}. Confirm by DNA sequencing.")
            recs.append("Confirm at DNA level (Sanger/NGS) to establish nucleotide change and zygosity.")
        elif "AR" in inh:
            if hom:
                status = "Biallelic (apparent homozygous) P/LP variant"
                lines.append(f"Apparent homozygous P/LP variant(s) in {gene} ({', '.join(v['hgvs_p'] or v['hgvs_c'] for v in hom)}): "
                             f"consistent with a diagnosis of {dnames} if the phenotype is concordant. "
                             "Exclude allele drop-out / a large deletion on the other allele (parental testing).")
            elif len(het) >= 2:
                status = "Two heterozygous P/LP variants (possible compound heterozygote)"
                lines.append(f"Two or more heterozygous P/LP variants in {gene}: possible compound heterozygosity. "
                             "Phase must be established (parental testing or long-read) before diagnosis.")
                recs.append("Parental testing to establish phase (cis vs trans).")
            else:
                status = "Heterozygous carrier of a P/LP variant (AR condition)"
                lines.append(f"Single heterozygous P/LP variant in {gene} associated with autosomal recessive "
                             f"{dnames}: consistent with carrier status. A second variant (including CNV or deep "
                             "intronic) should be sought if the phenotype is strongly suggestive.")
                recs.append("Offer partner / reproductive carrier testing.")
        elif inh.startswith("X"):
            status = "P/LP variant in X-linked gene"
            if sex in ("male", "m"):
                lines.append(f"Hemizygous P/LP variant in {gene} in a male: consistent with {dnames}.")
            else:
                lines.append(f"P/LP variant in X-linked {gene}: females may be carriers or variably affected "
                             "(X-inactivation); interpret with phenotype.")
            recs.append("Carrier testing of at-risk female relatives.")
        else:
            status = "P/LP variant identified"
            lines.append(f"P/LP variant(s) in {gene}: consistent with {dnames} (or predisposition) "
                         f"in the appropriate clinical context.")
        recs += ["Confirm P/LP variant(s) with an orthogonal method (e.g. Sanger/NGS on a new sample).",
                 "Genetic counselling; offer cascade testing to at-risk relatives."]
    else:
        lines.append(f"No pathogenic or likely pathogenic variants identified in the analysed region of {gene}.")
    if vus:
        lines.append(f"{len(vus)} variant(s) of uncertain significance - should not be used for clinical "
                     "decision-making; consider segregation/functional data and periodic reanalysis.")
        recs.append("Re-evaluate VUS periodically (e.g. every 1-2 years) or when new evidence emerges.")
    if ctx.get("acmg_sf") and plp:
        lines.append(f"{gene} is on the ACMG secondary-findings list ({ACMG_SF_VERSION}); P/LP findings are "
                     "medically actionable and reportable as secondary findings (per patient consent).")
    if sample["qc"].get("warnings"):
        recs.append("Review QC warnings before reporting.")
    recs += list(ctx.get("management") or [])
    recs = list(dict.fromkeys(recs))
    worst = None
    for c in CLASS_ORDER:
        if any(v["classification"] == c for v in vs):
            worst = c
            break
    return {"gene": gene, "inheritance": inh or "not specified", "diseases": diseases, "status": status,
            "highest_class": worst or "None", "interpretation": lines, "recommendations": recs,
            "acmg_sf_gene": bool(ctx.get("acmg_sf"))}


def most_severe(conss):
    for c in CONSEQUENCE_SEVERITY:
        if c in conss:
            return c
    return None
