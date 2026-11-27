"""Annotation loading and lookup (gene/transcript context, known variants, population data).

Annotation JSON schema (all keys optional; coordinates are 1-based, inclusive):

{
  "gene": "HBB", "transcript": "NM_000518.5", "protein_id": "NP_000509.1",
  "chrom": "chr11", "genomic_offset": 0,          # for VCF export
  "cds": [1, 444],                                  # CDS on the reference incl. stop codon
  "exons": [[1, 142], [273, 495]],                  # genomic-style references only
  "exon_boundaries_c": [92, 315],                   # cDNA references: last c. base of each exon
  "inheritance": "AR", "lof_mechanism": true,
  "missense_constrained": false, "truncating_only": false,
  "domains": [{"name": "Globin", "start": 3, "end": 147, "critical": false}],   # protein coords
  "hotspots": [7],                                  # protein residues (PM1)
  "known_variants": [{"hgvs_c": "c.20A>T", "hgvs_p": "p.Glu7Val", "classification": "Pathogenic",
                      "source": "ClinVar VCV000015333"}],
  "population": {"c.20A>T": {"af": 0.004, "source": "gnomAD"}},  # or {"c.20A>T": 0.004}
  "population_dataset": "gnomAD v4.1",              # declares population data was checked (enables PM2)
  "insilico": {"c.20A>T": {"REVEL": 0.85, "SpliceAI": 0.01}},
  "ba1_exceptions": ["c.20A>T"],
  "thresholds": {"ba1": 0.05, "bs1": 0.01, "pm2": 0.0001},
  "diseases": [...], "management": [...], "phenotype_keywords": [...],
  "references": {"REF_ID": { ...same keys, per reference sequence... }}
}
"""

from __future__ import annotations

import json
import re

from .genetic_code import AA1
from .knowledge import ACMG_SF_GENES, gene_info


def load_annotation(path):
    if not path:
        return {}
    with open(path) as fh:
        return json.load(fh)


def _merge(base, over):
    out = dict(base)
    for k, v in (over or {}).items():
        if k == "references":
            continue
        out[k] = v
    return out


def user_context(ann, ref_id):
    ann = ann or {}
    top = {k: v for k, v in ann.items() if k != "references"}
    refs = ann.get("references") or {}
    per = refs.get(ref_id)
    if per is None:  # tolerate missing/different version numbers (NM_000518 vs NM_000518.5)
        base = str(ref_id).split(".")[0]
        per = next((v for k, v in refs.items() if str(k).split(".")[0] == base), {})
    return _merge(top, per)


def _ref_keyed(ann, ref_id):
    refs = (ann or {}).get("references") or {}
    base = str(ref_id).split(".")[0]
    return any(str(k) == str(ref_id) or str(k).split(".")[0] == base for k in refs)


def context_for_reference(ann, ref_id, gene=None, auto=None, manual_known=None):
    """Precedence: online auto-annotation < built-in knowledge < user annotation file;
    known variants from all sources (plus manual entries) are combined."""
    ctx = user_context(ann, ref_id)
    symbol = gene or ctx.get("gene") or (auto or {}).get("gene")
    if symbol:
        # per-reference settings may also be keyed by gene symbol (--ref-setting CELF2:inheritance=AR);
        # settings keyed by the reference ID take precedence
        refs = (ann or {}).get("references") or {}
        by_gene = next((v for k, v in refs.items() if str(k).upper() == str(symbol).upper()), None)
        if by_gene and not _ref_keyed(ann, ref_id):
            ctx = _merge(ctx, by_gene)
    base = _merge(auto or {}, gene_info(symbol) if symbol else {})
    merged = _merge(base, ctx)
    kv = list((auto or {}).get("known_variants", [])) + list(ctx.get("known_variants", [])) + list(manual_known or [])
    seen, uniq = set(), []
    for x in reversed(kv):  # later (user / manual) entries win
        key = normalize_hgvs(x.get("hgvs_c", "")) + "|" + normalize_hgvs(x.get("hgvs_p", ""))
        if key not in seen:
            seen.add(key)
            uniq.append(x)
    merged["known_variants"] = list(reversed(uniq))
    if auto and auto.get("hotspots") and ctx.get("hotspots"):
        merged["hotspots"] = sorted(set(auto["hotspots"]) | set(ctx["hotspots"]))
    if symbol:
        merged["gene"] = symbol.upper()
    merged["acmg_sf"] = bool(merged.get("gene") and merged["gene"] in ACMG_SF_GENES)
    # normalise population and in-silico tables
    pop = {}
    for k, v in (merged.get("population") or {}).items():
        if isinstance(v, dict):
            pop[normalize_hgvs(k)] = v
        else:
            pop[normalize_hgvs(k)] = {"af": float(v)}
    merged["population"] = pop
    merged["insilico"] = {normalize_hgvs(k): v for k, v in (merged.get("insilico") or {}).items()}
    merged["ba1_exceptions"] = {normalize_hgvs(x) for x in merged.get("ba1_exceptions", [])}
    return merged


def normalize_hgvs(s: str) -> str:
    """Canonical form for matching: strip accession, parentheses, spaces, trailing del/dup bases."""
    if not s:
        return ""
    s = s.strip()
    if ":" in s:
        s = s.split(":", 1)[1]
    s = s.replace("(", "").replace(")", "").replace(" ", "")
    s = re.sub(r"del[ACGTN]+ins", "delins", s)  # c.470_471delGAinsAC == c.470_471delinsAC
    s = re.sub(r"(del|dup)[ACGTN]+$", r"\1", s)
    s = s.replace("*", "Ter") if s.startswith("p.") else s
    if s.startswith("p.") and re.match(r"^p\.[A-Z]\d+[A-Z]$", s):  # one-letter missense
        m = re.match(r"^p\.([A-Z])(\d+)([A-Z])$", s)
        from .genetic_code import AA3
        s = f"p.{AA3.get(m.group(1), 'Xaa')}{m.group(2)}{AA3.get(m.group(3), 'Xaa')}"
    return s


_P_MISSENSE = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)$")


def parse_missense(hgvs_p):
    """Return (pos, ref1, alt1) for a missense/nonsense p. string, else None."""
    m = _P_MISSENSE.match(normalize_hgvs(hgvs_p))
    if not m:
        return None
    r, pos, a = m.groups()
    return int(pos), AA1.get(r, "X"), AA1.get(a, "X")


def is_pathogenic_class(c):
    c = (c or "").lower().replace("_", " ")
    return c.startswith("pathogenic") or c.startswith("likely pathogenic") or c in ("p", "lp", "p/lp")


def population_lookup(ctx, keys):
    for k in keys:
        nk = normalize_hgvs(k)
        if nk and nk in ctx["population"]:
            return ctx["population"][nk]
    return None


def insilico_lookup(ctx, keys):
    for k in keys:
        nk = normalize_hgvs(k)
        if nk and nk in ctx["insilico"]:
            return ctx["insilico"][nk]
    return None
