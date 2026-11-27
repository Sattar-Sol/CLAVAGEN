"""Benchmark against the ClinGen Evidence Repository (expert-panel classifications).

For each expert-classified variant with a RefSeq transcript (NM_) c. description
in an exon (substitutions and small indels), a patient sequence carrying the variant
is built from the RefSeq transcript and analysed by CLAVAGEN exactly as a user
sequence would be. The ClinVar-based reported class is switched OFF, and the
variant's own record is removed from the known-variant list (leave-one-out), so the
result measures the automated ACMG engine only.

Outputs: per-variant TSV, summary JSON and an HTML report with agreement,
confusion matrix, critical (P/LP <-> B/LB) discordance and per-criterion
precision / recall against the expert panels' applied criteria.

Data source: ClinGen Evidence Repository (https://erepo.clinicalgenome.org), tabbed export.
"""

from __future__ import annotations

import csv
import datetime as _dt
import io
import json
import os
import random
import re
from html import escape

from . import __version__

EREPO_URL = "https://erepo.clinicalgenome.org/redmine/projects/evrepo/pcer/api/classifications/all?format=tabbed"
CLASSES = ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]
_POS = r"(?:-\d+|\*\d+|\d+)"
_HGVS = re.compile(rf"(N[MR]_\d+\.\d+)(?:\([A-Za-z0-9\-]+\))?:(c\.{_POS}(?:_{_POS})?[A-Za-z>0-9]*)")


def norm_class(s):
    s = (s or "").strip().lower()
    for c in CLASSES:
        if s == c.lower():
            return c
    if s.startswith("likely path"):
        return "Likely pathogenic"
    if s.startswith("path"):
        return "Pathogenic"
    if s.startswith("likely ben"):
        return "Likely benign"
    if s.startswith("ben"):
        return "Benign"
    if s.startswith("uncertain") or s == "vus":
        return "Uncertain significance"
    return ""


def tier(c):
    return {"Pathogenic": "P/LP", "Likely pathogenic": "P/LP", "Uncertain significance": "VUS",
            "Likely benign": "B/LB", "Benign": "B/LB"}.get(c, "")


def _col(header, *keys):
    for h in header:
        hl = h.lower().lstrip("#").strip()
        if all(k in hl for k in keys):
            return h
    return None


def parse_codes(s):
    """'PS3_Moderate, PM2_Supporting; PP1' -> {'PS3', 'PM2', 'PP1'} (base codes)."""
    out = set()
    for tok in re.split(r"[,;\s]+", s or ""):
        m = re.match(r"^(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])", tok.strip().upper())
        if m:
            out.add(m.group(1))
    return out


_STRENGTH = {"supporting": "Supporting", "moderate": "Moderate", "strong": "Strong", "verystrong": "VeryStrong",
             "standalone": "StandAlone"}

# Criteria a sequence-based tool can assess from the variant and public databases alone. The others need
# literature, functional, case, segregation or phenotype data (PS2-PS4, PM3, PM6, PP1, PP4, BS2-BS4, BP2, BP5).
COMPUTABLE = {"PVS1", "PS1", "PM1", "PM2", "PM4", "PM5", "PP2", "PP3", "BA1", "BS1", "BP1", "BP3", "BP4", "BP7"}


def parse_criteria(s):
    """'PS3_Moderate, PM2_Supporting; PP1' -> {'PS3': 'Moderate', 'PM2': 'Supporting', 'PP1': 'Supporting'}.
    A code without a strength gets its default ACMG/AMP strength."""
    from .acmg import CRITERIA
    out = {}
    for m in re.finditer(r"(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(?:[_ ]+(very[ _-]?strong|stand[ _-]?alone|"
                         r"strong|moderate|supporting))?", (s or "").upper(), re.I):
        code = m.group(1)
        st = m.group(2)
        st = _STRENGTH.get(re.sub(r"[ _-]", "", st.lower())) if st else None
        out[code] = st or CRITERIA.get(code, ("", "Supporting", ""))[1]
    return out


def combine_codes(criteria):
    """ACMG/AMP class from a {code: strength} dict, using the same combining rules as CLAVAGEN."""
    from .acmg import Criterion, combine
    return combine([Criterion(c, "met", st, "") for c, st in criteria.items()])["classification"]


def analyse_disagreement(row):
    """Why does CLAVAGEN disagree with the expert panel? Adds to ``row``:

    expert_computable_class - the expert classification recomputed from the computable criteria only
    disagreement            - '' (agree), 'missing evidence', 'criterion assignment' or 'evidence strength'
    extra_codes / missed_codes / strength_diffs / expert_noncomputable - the details
    """
    e = row.get("expert_criteria") or {c: None for c in row.get("expert_codes", [])}
    t = row.get("tool_criteria") or {}
    e_comp = {c: st for c, st in e.items() if c in COMPUTABLE}
    t_comp = {c: st for c, st in t.items() if c in COMPUTABLE}
    row["expert_noncomputable"] = sorted(c for c in e if c not in COMPUTABLE)
    row["extra_codes"] = sorted(set(t_comp) - set(e_comp))
    row["missed_codes"] = sorted(set(e_comp) - set(t_comp))
    row["strength_diffs"] = sorted(f"{c}: tool {t_comp[c]} / expert {e_comp[c]}" for c in set(t_comp) & set(e_comp)
                                   if e_comp[c] and t_comp[c] and e_comp[c] != t_comp[c])
    row["expert_computable_class"] = combine_codes({c: st or "Supporting" for c, st in e_comp.items()})
    # simulated: the user enters the clinical evidence the experts had (criteria a sequence tool cannot compute),
    # CLAVAGEN keeps its own automated criteria
    from .acmg import CRITERIA
    sup = dict(t)
    for c in row["expert_noncomputable"]:
        if c in CRITERIA:
            sup[c] = e.get(c) or CRITERIA[c][1]
    row["supplemented_class"] = combine_codes({c: st or "Supporting" for c, st in sup.items()})
    if row.get("tool_class") == row["expert_class"]:
        row["disagreement"] = ""
    elif row.get("tool_class") == row["expert_computable_class"]:
        # with the same kind of evidence the experts reach the tool's class: the gap is evidence a
        # sequence-based tool cannot have (functional, case, segregation, phenotype data)
        row["disagreement"] = "missing evidence"
    elif row["extra_codes"] or row["missed_codes"]:
        row["disagreement"] = "criterion assignment"
    elif row["strength_diffs"]:
        row["disagreement"] = "evidence strength"
    else:
        row["disagreement"] = "missing evidence"
    return row


def load_erepo(text):
    """Parse the Evidence Repository tabbed export into a list of usable records."""
    rows = list(csv.DictReader(io.StringIO(text.lstrip("﻿")), delimiter="\t"))
    if not rows:
        return [], {"total": 0}
    header = list(rows[0].keys())
    c_title = _col(header, "variation") or header[0]
    c_hgvs = _col(header, "hgvs")
    c_gene = _col(header, "gene")
    c_assert = _col(header, "assertion") or _col(header, "classification")
    c_met = next((h for h in header if "evidence" in h.lower() and "not met" not in h.lower() and "met" in h.lower()), None)
    c_moi = _col(header, "inheritance")
    c_panel = _col(header, "expert panel") or _col(header, "panel")
    c_retr = _col(header, "retract")
    c_car = _col(header, "allele registry")
    c_cv = _col(header, "clinvar")
    c_dis = _col(header, "disease")
    stats = {"total": len(rows), "retracted": 0, "no_class": 0, "no_transcript_hgvs": 0, "intronic_or_complex": 0,
             "used": 0}
    out = []
    for r in rows:
        if c_retr and str(r.get(c_retr, "")).strip().lower() in ("true", "yes", "1", "retracted"):
            stats["retracted"] += 1
            continue
        cls = norm_class(r.get(c_assert, ""))
        if not cls:
            stats["no_class"] += 1
            continue
        texts = [r.get(c_title, "") or "", r.get(c_hgvs, "") or ""]
        hits = [m for t in texts for m in _HGVS.finditer(t)]
        hits = [m for m in hits if m.group(1).startswith("NM_")]
        if not hits:
            stats["no_transcript_hgvs"] += 1
            continue
        tx, c = hits[0].group(1), hits[0].group(2)
        if re.search(r"\d[+-]\d", c) or not parse_c(c):
            stats["intronic_or_complex"] += 1
            continue
        gene = (r.get(c_gene) or "").strip()
        if not gene:
            m = re.search(r"\(([A-Za-z0-9\-]+)\):", r.get(c_title, "") or "")
            gene = m.group(1) if m else ""
        out.append({"title": (r.get(c_title) or "").strip(), "transcript": tx, "hgvs_c": c, "gene": gene,
                    "expert_class": cls, "expert_codes": sorted(parse_codes(r.get(c_met, "") if c_met else "")),
                    "expert_criteria": parse_criteria(r.get(c_met, "") if c_met else ""),
                    "moi": (r.get(c_moi) or "").strip() if c_moi else "",
                    "panel": (r.get(c_panel) or "").strip() if c_panel else "",
                    "caid": (r.get(c_car) or "").strip() if c_car else "",
                    "clinvar_id": (r.get(c_cv) or "").strip() if c_cv else "",
                    "disease": (r.get(c_dis) or "").strip() if c_dis else ""})
    stats["used"] = len(out)
    return out, stats


# ------------------------------------------------------------------ HGVS -> sequence
_C = re.compile(rf"^c\.({_POS})(?:_({_POS}))?(?:([ACGT])>([ACGT])|del([ACGT]*)ins([ACGT]+)|delins([ACGT]+)|"
                rf"del([ACGT]*)|dup([ACGT]*)|ins([ACGT]+))$")


def parse_c(c):
    m = _C.match(c)
    if not m:
        return None
    a, b = m.group(1), m.group(2) or m.group(1)
    if m.group(3):
        return ("sub", a, b, m.group(3), m.group(4))
    if m.group(6) or m.group(7):
        return ("delins", a, b, "", m.group(6) or m.group(7))
    if m.group(8) is not None and "del" in c and "ins" not in c:
        return ("del", a, b, m.group(8), "")
    if m.group(9) is not None and "dup" in c:
        return ("dup", a, b, m.group(9), "")
    if m.group(10):
        return ("ins", a, b, "", m.group(10))
    return None


def c_to_index(pos, cds_start, cds_end):
    """c. position -> 0-based index in the mRNA (cds_start/end 1-based inclusive)."""
    if pos.startswith("*"):
        return cds_end - 1 + int(pos[1:])
    p = int(pos)
    # c.1 = A of ATG (index cds_start-1); c.-1 = base before it (index cds_start-2)
    return cds_start - 2 + p if p > 0 else cds_start - 1 + p


def apply_c(mrna, cds, c):
    """Apply a c. variant to the mRNA. Returns (mutated, start0, end0) or raises ValueError."""
    pv = parse_c(c)
    if not pv:
        raise ValueError(f"unsupported HGVS {c}")
    kind, a, b, ref, alt = pv
    cs, ce = cds
    i, j = c_to_index(a, cs, ce), c_to_index(b, cs, ce)
    if not (0 <= i < len(mrna) and 0 <= j < len(mrna)):
        raise ValueError(f"{c} outside transcript")
    if kind == "sub":
        if mrna[i] != ref:
            raise ValueError(f"reference base mismatch for {c} (transcript has {mrna[i]})")
        return mrna[:i] + alt + mrna[i + 1:], i, i + 1
    if kind == "del":
        return mrna[:i] + mrna[j + 1:], i, j + 1
    if kind == "dup":
        return mrna[:j + 1] + mrna[i:j + 1] + mrna[j + 1:], i, j + 1
    if kind == "ins":
        return mrna[:i + 1] + alt + mrna[j:], i, j
    if kind == "delins":
        return mrna[:i] + alt + mrna[j + 1:], i, j + 1
    raise ValueError(c)


# ------------------------------------------------------------------ run
VARIANT_TYPE = {"missense_variant": "missense", "synonymous_variant": "synonymous", "stop_gained": "nonsense",
                "frameshift_variant": "frameshift", "inframe_deletion": "inframe indel",
                "inframe_insertion": "inframe indel", "start_lost": "start/stop loss", "stop_lost": "start/stop loss",
                "splice_region_variant": "splice region", "5_prime_UTR_variant": "UTR", "3_prime_UTR_variant": "UTR"}


NETWORK_ERR = r"(no internet connection|timed out|HTTP 5\d\d|HTTP 429|Connection reset|Temporary failure)"


def _key(r):
    return f"{r['transcript']}:{r['hgvs_c']}"


def run_benchmark(records, out_dir, flank=150, max_n=None, seed=7, genes=None, log=None, online=True,
                  resume=False, meta=None, figures=True):
    """Analyse each record; per-variant results are appended to benchmark_results.jsonl as they finish,
    so an interrupted run continues with ``resume=True`` (same records, seed and --max)."""
    import time
    from .annotation import normalize_hgvs
    from .online import fetch_refseq
    from .pipeline import Pipeline
    log = log or print
    os.makedirs(out_dir, exist_ok=True)
    recs = [r for r in records if not genes or r["gene"].upper() in {g.upper() for g in genes}]
    if max_n and len(recs) > max_n:
        rnd = random.Random(seed)
        by = {}
        for r in recs:
            by.setdefault(tier(r["expert_class"]), []).append(r)
        per = max(1, max_n // max(1, len(by)))
        pick = []
        for t, lst in by.items():
            rnd.shuffle(lst)
            pick += lst[:per]
        recs = pick[:max_n]
        log(f"Stratified sample of {len(recs)} variants (seed {seed}).")
    tx_cache, auto_cache, results = {}, {}, []
    ref_dir = os.path.join(out_dir, "references")
    os.makedirs(ref_dir, exist_ok=True)
    jl = os.path.join(out_dir, "benchmark_results.jsonl")
    done = {}
    if resume and os.path.isfile(jl):
        with open(jl) as fh:
            for line in fh:
                try:
                    d = json.loads(line)
                    if d.get("status") != "error":      # errors (e.g. internet drop) are tried again
                        done[_key(d)] = d
                except (ValueError, KeyError):
                    pass  # a line cut off by the interruption
        log(f"Resuming: {len(done)} variant(s) already analysed in {jl}.")
    meta = dict(meta or {})
    meta.update({"started": meta.get("started") or _dt.datetime.now().isoformat(timespec="seconds"),
                 "flank": flank, "max": max_n or 0, "seed": seed, "genes": genes or "all", "online": online,
                 "selected": len(recs), "resumed_from": len(done)})
    t_start = time.time()
    jfh = open(jl, "a" if resume else "w")
    for k, r in enumerate(recs, 1):
        if _key(r) in done:
            results.append(done[_key(r)])
            continue
        t0 = time.time()
        row = dict(r, tool_class="", tool_acmg="", tool_codes=[], status="", note="")
        try:
            tx = r["transcript"]
            if tx not in tx_cache:
                tx_cache[tx] = fetch_refseq(tx)
            gb = tx_cache[tx]
            if not gb.get("cds"):
                raise ValueError("no CDS in RefSeq record")
            mut, s0, e0 = apply_c(gb["seq"], gb["cds"], r["hgvs_c"])
            lo, hi = max(0, s0 - flank), min(len(mut), s0 + flank + max(0, len(mut) - len(gb["seq"])) + 1)
            query = mut[lo:hi]
            ref_path = os.path.join(ref_dir, f"{tx}.fa")
            if not os.path.isfile(ref_path):
                with open(ref_path, "w") as fh:
                    fh.write(f">{tx} ({gb.get('gene') or r['gene']}) RefSeq\n{gb['seq']}\n")
            pl = Pipeline([f"bench{k}={query}"], [ref_path], gene=gb.get("gene") or r["gene"], online=online,
                          use_clinvar=False, auto_cache=auto_cache, log=lambda *_: None,
                          exclude_known={normalize_hgvs(r["hgvs_c"])} | (
                              {normalize_hgvs(m.group(1))} if (m := re.search(r"\((p\.[^)]*)\)\s*$", r.get("title", "")))
                              else set()))
            res = pl.run()
            vs = res["samples"][0]["variants"]
            target = [v for v in vs if normalize_hgvs(v.get("hgvs_c", "")) == normalize_hgvs(r["hgvs_c"])]
            if not target:
                # same change written differently (e.g. c.904delT vs 3'-shifted c.907del): compare the sequences
                for v in vs:
                    try:
                        if v.get("hgvs_c") and apply_c(gb["seq"], gb["cds"], v["hgvs_c"])[0] == mut:
                            target = [v]
                            row["note"] = f"expert HGVS {r['hgvs_c']} = {v['hgvs_c']} (HGVS 3' rule)"
                            break
                    except Exception:
                        pass
            if not target:
                row["status"] = "variant not recovered"
                row["note"] = "called: " + ", ".join(v.get("hgvs_c", "") for v in vs)
            else:
                v = target[0]
                row.update(tool_class=v["classification"], tool_acmg=v.get("acmg_classification", v["classification"]),
                           tool_codes=sorted({c.split("_")[0] for c in v["criteria_met"]}), status="ok",
                           tool_criteria=parse_criteria(" ".join(v["criteria_met"])), points=v.get("points"))
                row["population_af"] = v.get("population_af")
                row["revel"] = (v.get("insilico") or {}).get("REVEL")
                row["spliceai"] = (v.get("insilico") or {}).get("SpliceAI")
                row["hgvs_p"] = v.get("hgvs_p", "")
                row["consequence"] = v.get("consequence", "")
                row["variant_type"] = VARIANT_TYPE.get(v.get("consequence", ""), v.get("consequence", "") or "other")
                row["tool_criteria_str"] = ",".join(f"{c}_{st}" for c, st in sorted(row["tool_criteria"].items()))
                analyse_disagreement(row)
                if len(vs) > 1:
                    row["note"] = "; ".join(filter(None, [row["note"], f"{len(vs) - 1} additional call(s) in the window"]))
            # keep data-source problems visible (e.g. old transcript version, VEP lookup failure)
            qn = [n for n in (res["samples"][0].get("qc") or {}).get("notes", [])
                  if "VEP" in n or "RefSeq version" in n or "current version" in n]
            if qn:
                row["note"] = "; ".join(filter(None, [row["note"]] + qn))[:400]
            if row["status"] == "ok":
                row["lookup_transcript"] = v.get("lookup_transcript", "")
                row["vep_resolved"] = v.get("vep_resolved")
            allnotes = " ".join((res["samples"][0].get("qc") or {}).get("notes", []))
            if row["status"] == "ok" and re.search(NETWORK_ERR, allnotes):
                # a dropped connection is not a property of the variant: mark it so --resume tries again
                row["status"] = "error"
                row["note"] = "network: " + re.search(NETWORK_ERR + r"[^;]*", allnotes).group(0)[:160]
        except Exception as e:
            row["status"] = "error"
            row["note"] = str(e)[:200]
        row["seconds"] = round(time.time() - t0, 2)
        results.append(row)
        jfh.write(json.dumps(row, default=str) + "\n")
        jfh.flush()
        if k % 10 == 0 or k == len(recs):
            el = time.time() - t_start
            todo = len(recs) - k
            new = max(1, k - len([1 for x in recs[:k] if _key(x) in done]))
            eta = f", about {el / new * todo / 60:.0f} min left" if todo and new >= 10 else ""
            log(f"  {k}/{len(recs)} variants processed{eta}")
    jfh.close()
    meta["finished"] = _dt.datetime.now().isoformat(timespec="seconds")
    meta["run_seconds"] = round(time.time() - t_start)
    summary = summarise(results, meta=meta, seed=seed)
    write_outputs(results, summary, out_dir, figures=figures, log=log)
    return summary


NOT_TESTABLE = r"(outside transcript|reference base mismatch)"


def _not_testable(r):
    """The expert HGVS cannot be placed on the RefSeq transcript record (e.g. a promoter position beyond
    its 5' end, or a reference base that differs): a property of the input, not a tool error."""
    return r["status"] == "error" and bool(re.search(NOT_TESTABLE, r.get("note", "")))


def summarise(results, meta=None, seed=7):
    from .benchmark_stats import publication_stats
    ok = [r for r in results if r["status"] == "ok"]
    n = len(ok)
    exact = sum(1 for r in ok if r["tool_class"] == r["expert_class"])
    t3 = sum(1 for r in ok if tier(r["tool_class"]) == tier(r["expert_class"]))
    crit = sum(1 for r in ok if {tier(r["tool_class"]), tier(r["expert_class"])} == {"P/LP", "B/LB"})
    conf = {e: {t: 0 for t in CLASSES} for e in CLASSES}
    for r in ok:
        conf[r["expert_class"]][r["tool_class"]] += 1
    sens = {}
    for t in ("P/LP", "VUS", "B/LB"):
        tot = [r for r in ok if tier(r["expert_class"]) == t]
        sens[t] = {"n": len(tot), "agree": sum(1 for r in tot if tier(r["tool_class"]) == t)}
    codes = {}
    for r in ok:
        e, m = set(r["expert_codes"]), set(r["tool_codes"])
        for c in e | m:
            d = codes.setdefault(c, {"tp": 0, "fp": 0, "fn": 0})
            if c in e and c in m:
                d["tp"] += 1
            elif c in m:
                d["fp"] += 1
            else:
                d["fn"] += 1
    for c, d in codes.items():
        d["precision"] = round(d["tp"] / (d["tp"] + d["fp"]), 3) if d["tp"] + d["fp"] else None
        d["recall"] = round(d["tp"] / (d["tp"] + d["fn"]), 3) if d["tp"] + d["fn"] else None
    kinds = ("missing evidence", "criterion assignment", "evidence strength")
    disc = [r for r in ok if r["tool_class"] != r["expert_class"]]
    for r in disc:
        if "disagreement" not in r:
            analyse_disagreement(r)
    disagreement = {k: sum(1 for r in disc if r.get("disagreement") == k) for k in kinds}
    for r in ok:
        if "expert_computable_class" not in r:
            analyse_disagreement(r)
    comp_exact = sum(1 for r in ok if r["tool_class"] == r["expert_computable_class"])
    comp_t3 = sum(1 for r in ok if tier(r["tool_class"]) == tier(r["expert_computable_class"]))
    genes = {}
    for r in ok:
        g = genes.setdefault(r["gene"], {"n": 0, "agree3": 0})
        g["n"] += 1
        g["agree3"] += tier(r["tool_class"]) == tier(r["expert_class"])
    return {"software": f"CLAVAGEN {__version__}", "date": _dt.datetime.now().isoformat(timespec="seconds"),
            "attempted": len(results), "analysed": n,
            "not_recovered": sum(1 for r in results if r["status"] == "variant not recovered"),
            "not_testable": sum(1 for r in results if _not_testable(r)),
            "errors": sum(1 for r in results if r["status"] == "error" and not _not_testable(r)),
            "exact_agreement": round(exact / n, 3) if n else None,
            "three_tier_agreement": round(t3 / n, 3) if n else None,
            "critical_discordance": round(crit / n, 3) if n else None, "critical_n": crit,
            "tier_agreement": sens, "confusion": conf, "criteria": dict(sorted(codes.items())), "genes": genes,
            "discordant": len(disc), "disagreement": disagreement,
            "computable_exact_agreement": round(comp_exact / n, 3) if n else None,
            "computable_three_tier_agreement": round(comp_t3 / n, 3) if n else None,
            "computable_criteria": sorted(COMPUTABLE),
            "statistics": publication_stats(ok, seed=seed) if ok else {},
            "status_counts": {s: sum(1 for r in results if r["status"] == s) for s in sorted({r["status"] for r in results})},
            "run": meta or {}}


TSV_COLUMNS = ["gene", "transcript", "hgvs_c", "hgvs_p", "variant_type", "expert_class", "tool_class", "agree_5class",
               "agree_3tier", "expert_codes", "tool_codes", "expert_criteria", "tool_criteria", "points",
               "expert_computable_class", "disagreement", "tool_only_codes", "expert_only_codes", "strength_differences",
               "expert_noncomputable_codes", "population_af", "revel", "spliceai", "vep_resolved", "lookup_transcript", "panel", "disease",
               "moi", "caid", "clinvar_id", "seconds", "status", "note"]


def _tsv_row(r):
    ok = r["status"] == "ok"
    ec = r.get("expert_criteria") or {}
    val = {"agree_5class": r.get("tool_class") == r["expert_class"] if ok else "",
           "agree_3tier": tier(r.get("tool_class", "")) == tier(r["expert_class"]) if ok else "",
           "expert_codes": ",".join(r.get("expert_codes", [])), "tool_codes": ",".join(r.get("tool_codes", [])),
           "expert_criteria": ",".join(f"{c}_{s}" if s else c for c, s in sorted(ec.items())),
           "tool_criteria": r.get("tool_criteria_str", ""),
           "tool_only_codes": ",".join(r.get("extra_codes", [])), "expert_only_codes": ",".join(r.get("missed_codes", [])),
           "strength_differences": "; ".join(r.get("strength_diffs", [])),
           "expert_noncomputable_codes": ",".join(r.get("expert_noncomputable", []))}
    return ["" if (x := val.get(c, r.get(c, ""))) is None else x for c in TSV_COLUMNS]


def _write_tables(summary, out_dir):
    """CSV tables ready for a manuscript (Table 1-4) and a methods/results text with the numbers filled in."""
    st = summary.get("statistics") or {}
    if not st:
        return
    td = os.path.join(out_dir, "tables")
    os.makedirs(td, exist_ok=True)

    def ci(d):
        return f"{100 * d['p']:.1f} ({100 * d['ci95'][0]:.1f}-{100 * d['ci95'][1]:.1f})" if d and d.get("p") is not None else ""

    def kci(d):
        return f"{d['kappa']:.3f} ({d['ci95'][0]:.3f}-{d['ci95'][1]:.3f})" if d and d.get("ci95", [None])[0] is not None else ""
    with open(os.path.join(td, "table1_overall.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Measure", "n / N", "% (95% CI)"])
        for lab, key in (("3-tier agreement (P/LP, VUS, B/LB)", "three_tier"), ("Exact 5-class agreement", "exact"),
                         ("Critical discordance (P/LP vs B/LB)", "critical_discordance"),
                         ("3-tier agreement, expert class from computable criteria", "computable_three_tier"),
                         ("Exact agreement, expert class from computable criteria", "computable_exact")):
            d = st[key]
            w.writerow([lab, f"{d['k']}/{d['n']}", ci(d)])
        for t in ("P/LP", "B/LB"):
            for m, lab in (("sensitivity", "sensitivity"), ("ppv", "positive predictive value")):
                d = st["detection"][t][m]
                w.writerow([f"{t} {lab}", f"{d['k']}/{d['n']}", ci(d)])
        w.writerow(["Cohen's kappa, 3-tier", "", kci(st["kappa_three_tier"])])
        w.writerow(["Linearly weighted kappa, 5-class", "", kci(st["weighted_kappa_five_class"])])
        sp = st.get("with_clinical_evidence") or {}
        if sp:
            for lab, key in (("With the experts' clinical evidence entered (simulation): 3-tier agreement", "three_tier"),
                             ("With the experts' clinical evidence entered (simulation): exact agreement", "exact"),
                             ("With the experts' clinical evidence entered (simulation): critical discordance",
                              "critical_discordance")):
                d = sp[key]
                w.writerow([lab, f"{d['k']}/{d['n']}", ci(d)])
            w.writerow(["With the experts' clinical evidence entered (simulation): Cohen's kappa, 3-tier", "",
                        kci(sp["kappa_three_tier"])])
    for name, key in (("table2_by_panel.csv", "by_panel"), ("table2b_by_gene.csv", "by_gene"),
                      ("table2c_by_variant_type.csv", "by_variant_type"), ("table2d_by_expert_tier.csv", "by_expert_tier")):
        with open(os.path.join(td, name), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Group", "n", "Expert P/LP", "Expert VUS", "Expert B/LB", "3-tier agreement % (95% CI)",
                        "Exact agreement % (95% CI)", "Critical discordance n"])
            for g, d in st[key].items():
                w.writerow([g, d["n"], *[d["expert_tiers"][t] for t in ("P/LP", "VUS", "B/LB")], ci(d["three_tier"]),
                            ci(d["exact"]), d["critical"]["k"]])
    with open(os.path.join(td, "table3_criteria.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Criterion", "Computable by CLAVAGEN", "Both (TP)", "Tool only (FP)", "Experts only (FN)",
                    "Precision % (95% CI)", "Recall % (95% CI)"])
        for c, d in summary["criteria"].items():
            w.writerow([c, "yes" if c in COMPUTABLE else "no", d["tp"], d["fp"], d["fn"],
                        ci(prop_(d["tp"], d["tp"] + d["fp"])), ci(prop_(d["tp"], d["tp"] + d["fn"]))])
    with open(os.path.join(td, "table4_disagreement.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Main cause of disagreement", "n", "% of discordant"])
        nd = summary.get("discordant") or 0
        for k, v in summary["disagreement"].items():
            w.writerow([k, v, f"{100 * v / nd:.1f}" if nd else ""])
    with open(os.path.join(td, "confusion_matrix.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Expert \\ CLAVAGEN", *CLASSES])
        for e in CLASSES:
            w.writerow([e, *[summary["confusion"][e][t] for t in CLASSES]])
    _write_methods(summary, out_dir, ci, kci)


def prop_(k, n):
    from .benchmark_stats import prop
    return prop(k, n)


def _write_methods(summary, out_dir, _ci, _kci):
    def ci(d):  # "66.7% (95% CI 52.5-78.3)"
        return (f"{100 * d['p']:.1f}% (95% CI {100 * d['ci95'][0]:.1f}-{100 * d['ci95'][1]:.1f})"
                if d and d.get("p") is not None else "n/a")

    def kci(d):
        k = _kci(d)
        return f"{k.split(' ')[0]} (95% CI {k.split('(')[1].rstrip(')')}), {_lk(d)}" if k else "not estimable (too few variants)"
    st, run = summary["statistics"], summary.get("run") or {}
    f = run.get("erepo_filter") or {}
    cov = st["evidence_coverage"]
    nc = st["noncomputable_criteria"]
    sample = (f", of which a stratified sample of {run.get('selected')} was analysed (seed {run.get('seed')})"
              if run.get("max") else "")
    gsel = "" if run.get("genes") in (None, "all") else f" Analysis was restricted to {', '.join(run['genes'])}."
    lines = [
        "# Benchmark: methods and results (numbers filled in by CLAVAGEN)", "",
        "Draft text for a manuscript. Check every sentence against your own understanding before using it.", "",
        "## Methods", "",
        f"Expert-panel classifications were taken from the ClinGen Evidence Repository export "
        f"({run.get('erepo_file', 'downloaded file')}, SHA-256 {run.get('erepo_sha256', 'n/a')[:16]}..., "
        f"{f.get('total', '?')} rows, retrieved {run.get('erepo_date', 'n/a')}). Retracted assertions "
        f"({f.get('retracted', '?')}), rows without a recognised classification ({f.get('no_class', '?')}), rows without an NM_ "
        f"transcript c. description ({f.get('no_transcript_hgvs', '?')}) and intronic or complex variants "
        f"({f.get('intronic_or_complex', '?')}) were excluded, leaving {f.get('used', '?')} exonic substitutions and small "
        f"indels{sample}.{gsel}",
        "",
        f"For each variant a patient sequence (variant +/- {run.get('flank', 150)} nt) was generated from the RefSeq "
        f"transcript and analysed with {summary['software']} exactly as user input: alignment to the transcript, "
        "variant calling, HGVS normalisation and automated ACMG/AMP classification with the ClinGen SVI refinements. "
        "Population frequencies (gnomAD v4 exomes and genomes) and in-silico scores (REVEL, SpliceAI) were obtained "
        "through the Ensembl VEP REST service; variants described on an older RefSeq version were looked up on the "
        "current version when the coding sequence length was unchanged. To avoid information leakage, the "
        "ClinVar-reported classification was not used, the variant's own ClinVar records were excluded "
        "(leave-one-out), and gene-level evidence derived from ClinVar (mutational hotspots, loss-of-function "
        "mechanism) was recounted without the variant.",
        "",
        "Agreement was measured on the 3-tier scale (P/LP, VUS, B/LB) and the 5-class scale. Proportions are given "
        "with Wilson score 95% confidence intervals; agreement beyond chance with Cohen's kappa (3-tier) and linearly "
        "weighted kappa (5-class), with percentile bootstrap 95% intervals "
        f"({st['kappa_three_tier'].get('bootstrap', 1000)} resamples, seed {st['kappa_three_tier'].get('seed', 7)}). "
        "Each disagreement was assigned one main cause: missing evidence (the expert class recalculated from the "
        "computable criteria alone equals the tool's class), criterion assignment, or evidence strength.",
        "", "## Results", "",
        f"Of {summary['attempted']} variants, {summary['analysed']} were analysed "
        f"({summary.get('not_testable', 0)} could not be placed on the RefSeq transcript record, "
        f"{summary['not_recovered']} were not recovered as a single call, {summary['errors']} failed). CLAVAGEN agreed with the expert "
        f"panels on the 3-tier scale for {ci(st['three_tier'])} of variants (Cohen's kappa {kci(st['kappa_three_tier'])}) "
        f"and on the exact 5-class scale for {ci(st['exact'])} (weighted kappa "
        f"{kci(st['weighted_kappa_five_class'])}). Critical discordance (P/LP versus B/LB) occurred in "
        f"{st['critical_discordance']['k']} variants ({ci(st['critical_discordance'])}).",
        "",
        f"Experts used at least one criterion that a sequence-based tool cannot evaluate (e.g. PP4, PM3, PS3, PS4) for "
        f"{ci(nc['variants_with_any'])} of variants. When the expert classification was recalculated from the "
        f"computable criteria only, 3-tier agreement was {ci(st['computable_three_tier'])}. Of "
        f"{summary['discordant']} discordant variants, " + ", ".join(
            f"{v} ({100 * v / summary['discordant']:.0f}%) were due to {k}" for k, v in summary["disagreement"].items()
            if summary["discordant"]) + ".",
        "",
        *([f"In a simulation of the guided clinical-evidence input, the criteria that the experts based on clinical, "
           f"functional, segregation or case data were added to CLAVAGEN's automated criteria at the experts' strengths. "
           f"3-tier agreement then rose to {ci(st['with_clinical_evidence']['three_tier'])} (Cohen's kappa "
           f"{kci(st['with_clinical_evidence']['kappa_three_tier'])}), with "
           f"{st['with_clinical_evidence']['critical_discordance']['k']} critical discordances. This is an upper bound: "
           f"it assumes the user has the same clinical evidence as the panel. The remaining disagreements reflect "
           f"gene-specific rules of the expert panels.", ""] if st.get("with_clinical_evidence") else []),
        f"Evidence coverage: gnomAD frequency for {ci(cov['population_af'])} of variants, REVEL for "
        f"{ci(cov['revel_missense'])} of missense variants; {cov['old_transcript_remapped']['k']} variants were looked up "
        f"on a newer transcript version.",
        "", f"Run: {run.get('started', '')} to {run.get('finished', '')} ({run.get('run_seconds', 0) // 60} min)"
        + (f", Python {run['python']} on {run.get('platform', '')}" if run.get("python") else "")
        + f"; Ensembl release {run.get('ensembl_release') or 'n/a'}.",
    ]
    with open(os.path.join(out_dir, "methods_and_results.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _lk(d):
    from .benchmark_stats import landis_koch
    return landis_koch(d.get("kappa")) + " agreement"


def _supp_html(summary, pct):
    sp = (summary.get("statistics") or {}).get("with_clinical_evidence") or {}
    if not sp:
        return ""
    k = sp["kappa_three_tier"].get("kappa")
    return (f"<h2>If the clinical evidence is entered</h2><p class='muted'>Simulation of the guided clinical-evidence "
            f"input: the criteria the experts based on clinical, functional, segregation or case data "
            f"({sp['variants_with_clinical_evidence']} variants) are added to CLAVAGEN's automated criteria at the experts' "
            f"strengths. An upper bound - it assumes the same clinical evidence as the panel.</p>"
            f"<div class='kpi'><div><b>{pct(sp['three_tier']['p'])}</b>3-tier agreement</div>"
            f"<div><b>{pct(sp['exact']['p'])}</b>exact 5-class agreement</div>"
            f"<div><b>{'' if k is None else f'{k:.2f}'}</b>Cohen's kappa (3-tier)</div>"
            f"<div><b>{sp['critical_discordance']['k']}</b>critical discordances</div></div>")


def write_outputs(results, summary, out_dir, figures=True, log=None):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "benchmark_results.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(TSV_COLUMNS)
        for r in results:
            w.writerow(_tsv_row(r))
    with open(os.path.join(out_dir, "benchmark_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    _write_tables(summary, out_dir)
    if figures and summary.get("statistics"):
        from .benchmark_figures import make_figures
        try:
            made = make_figures(summary, results, os.path.join(out_dir, "figures"))
            if log and made:
                log(f"Figures: {len(made)} written to {os.path.join(out_dir, 'figures')}")
        except Exception as e:  # figures must never cost the results of a long run
            if not isinstance(e, ImportError):
                if log:
                    log(f"Figures failed ({e}); results are saved. Re-draw with 'clavagen benchmark-figures {out_dir}'.")
                return
            if log:
                log("Figures skipped: matplotlib is not installed (pip install matplotlib, then "
                    "'clavagen benchmark-figures " + out_dir + "').")
    pct = lambda x: "-" if x is None else f"{100 * x:.1f}%"  # noqa: E731
    conf_rows = "".join(
        f"<tr><th>{escape(e)}</th>" + "".join(
            f"<td style='{'background:#e7f5ec;font-weight:bold' if e == t else ''}'>{summary['confusion'][e][t]}</td>"
            for t in CLASSES) + "</tr>" for e in CLASSES)
    crit_rows = "".join(
        f"<tr><td>{escape(c)}</td><td>{d['tp']}</td><td>{d['fp']}</td><td>{d['fn']}</td><td>{pct(d['precision'])}</td>"
        f"<td>{pct(d['recall'])}</td></tr>" for c, d in summary["criteria"].items())
    gene_rows = "".join(f"<tr><td>{escape(g)}</td><td>{d['n']}</td><td>{pct(d['agree3'] / d['n'] if d['n'] else None)}</td></tr>"
                        for g, d in sorted(summary["genes"].items(), key=lambda kv: -kv[1]["n"]))
    nd = summary.get("discordant") or 0
    dis_rows = "".join(f"<tr><td>{escape(k)}</td><td>{v}</td><td>{pct(v / nd if nd else None)}</td></tr>"
                       for k, v in (summary.get("disagreement") or {}).items())
    tiers = "".join(f"<tr><td>{t}</td><td>{d['n']}</td><td>{pct(d['agree'] / d['n'] if d['n'] else None)}</td></tr>"
                    for t, d in summary["tier_agreement"].items())
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>CLAVAGEN benchmark</title>
<style>body{{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:28px;color:#1d2433;max-width:1100px}}
table{{border-collapse:collapse;margin:8px 0 22px}}td,th{{border:1px solid #e5e7eb;padding:6px 10px;text-align:left}}th{{background:#f3f4f6}}
.kpi{{display:flex;gap:14px;flex-wrap:wrap}}.kpi div{{background:#f4f6f9;border-radius:10px;padding:12px 16px;min-width:150px}}
.kpi b{{display:block;font-size:24px;color:#13294b}}.muted{{color:#6b7280;font-size:13px}}</style></head><body>
<h1>Benchmark against the ClinGen Evidence Repository</h1>
<p class="muted">{escape(summary['software'])} | {escape(summary['date'])} | automated ACMG engine only (ClinVar-based reported class OFF,
the variant's own ClinVar record excluded, and ClinVar-derived hotspots / loss-of-function mechanism recounted without it) | patient sequences built from RefSeq transcripts</p>
<div class="kpi"><div><b>{summary['analysed']}</b>variants analysed</div><div><b>{pct(summary['three_tier_agreement'])}</b>3-tier agreement
(P/LP, VUS, B/LB)</div><div><b>{pct(summary['exact_agreement'])}</b>exact 5-class agreement</div>
<div><b>{pct(summary['critical_discordance'])}</b>critical discordance (P/LP vs B/LB; n={summary['critical_n']})</div></div>
<p class="muted">Attempted {summary['attempted']}; not testable on the RefSeq transcript {summary.get('not_testable', 0)}; not recovered {summary['not_recovered']}; errors {summary['errors']}.</p>
<h2>Why the classifications differ</h2>
<p class="muted">Each of the {summary.get('discordant', 0)} discordant variants is assigned one main cause.
<b>Missing evidence</b>: recomputing the expert classification from the computable criteria only gives the tool's class -
the experts used evidence a sequence-based tool cannot have (functional, case, segregation or phenotype data).
<b>Criterion assignment</b>: the tool applied a computable criterion the experts did not, or missed one they applied.
<b>Evidence strength</b>: same computable criteria, different strength.</p>
<table><tr><th>Cause</th><th>n</th><th>% of discordant</th></tr>{dis_rows}</table>
<div class="kpi"><div><b>{pct(summary.get('computable_three_tier_agreement'))}</b>3-tier agreement with the expert class
recomputed from computable criteria</div><div><b>{pct(summary.get('computable_exact_agreement'))}</b>exact 5-class agreement
(computable criteria)</div></div>
<p class="muted">Computable criteria: {escape(', '.join(summary.get('computable_criteria', [])))}.</p>
{_supp_html(summary, pct)}
<h2>Agreement by expert tier</h2><table><tr><th>Expert tier</th><th>n</th><th>Tool agrees</th></tr>{tiers}</table>
<h2>Confusion matrix</h2><p class="muted">Rows: expert panel; columns: CLAVAGEN.</p>
<table><tr><th></th>{''.join(f'<th>{escape(c)}</th>' for c in CLASSES)}</tr>{conf_rows}</table>
<h2>Criteria compared with expert panels</h2><p class="muted">TP: applied by both; FP: applied by the tool only; FN: applied by
the expert panel only. Literature-based criteria (PS3, PS4, PP1...) are expected to be missed by sequence-only automation.</p>
<table><tr><th>Criterion</th><th>TP</th><th>FP</th><th>FN</th><th>Precision</th><th>Recall</th></tr>{crit_rows}</table>
<h2>By gene</h2><table><tr><th>Gene</th><th>n</th><th>3-tier agreement</th></tr>{gene_rows}</table>
<p class="muted">Data: ClinGen Evidence Repository (expert-panel curated classifications). Per-variant results: benchmark_results.tsv.</p>
</body></html>"""
    with open(os.path.join(out_dir, "benchmark_report.html"), "w") as fh:
        fh.write(html)
