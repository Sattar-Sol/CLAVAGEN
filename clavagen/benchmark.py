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
    stats = {"total": len(rows), "retracted": 0, "no_transcript_hgvs": 0, "intronic_or_complex": 0, "used": 0}
    out = []
    for r in rows:
        if c_retr and str(r.get(c_retr, "")).strip().lower() in ("true", "yes", "1", "retracted"):
            stats["retracted"] += 1
            continue
        cls = norm_class(r.get(c_assert, ""))
        if not cls:
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
                    "moi": (r.get(c_moi) or "").strip() if c_moi else "",
                    "panel": (r.get(c_panel) or "").strip() if c_panel else ""})
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
def run_benchmark(records, out_dir, flank=150, max_n=None, seed=7, genes=None, log=None, online=True):
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
    for k, r in enumerate(recs, 1):
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
                          exclude_known={normalize_hgvs(r["hgvs_c"])})
            res = pl.run()
            vs = res["samples"][0]["variants"]
            target = [v for v in vs if normalize_hgvs(v.get("hgvs_c", "")) == normalize_hgvs(r["hgvs_c"])]
            if not target:
                row["status"] = "variant not recovered"
                row["note"] = "called: " + ", ".join(v.get("hgvs_c", "") for v in vs)
            else:
                v = target[0]
                row.update(tool_class=v["classification"], tool_acmg=v.get("acmg_classification", v["classification"]),
                           tool_codes=sorted({c.split("_")[0] for c in v["criteria_met"]}), status="ok",
                           points=v.get("points"))
                if len(vs) > 1:
                    row["note"] = f"{len(vs) - 1} additional call(s) in the window"
        except Exception as e:
            row["status"] = "error"
            row["note"] = str(e)[:200]
        results.append(row)
        if k % 10 == 0 or k == len(recs):
            log(f"  {k}/{len(recs)} variants processed")
    summary = summarise(results)
    write_outputs(results, summary, out_dir)
    return summary


def summarise(results):
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
    genes = {}
    for r in ok:
        g = genes.setdefault(r["gene"], {"n": 0, "agree3": 0})
        g["n"] += 1
        g["agree3"] += tier(r["tool_class"]) == tier(r["expert_class"])
    return {"software": f"CLAVAGEN {__version__}", "date": _dt.datetime.now().isoformat(timespec="seconds"),
            "attempted": len(results), "analysed": n,
            "not_recovered": sum(1 for r in results if r["status"] == "variant not recovered"),
            "errors": sum(1 for r in results if r["status"] == "error"),
            "exact_agreement": round(exact / n, 3) if n else None,
            "three_tier_agreement": round(t3 / n, 3) if n else None,
            "critical_discordance": round(crit / n, 3) if n else None, "critical_n": crit,
            "tier_agreement": sens, "confusion": conf, "criteria": dict(sorted(codes.items())), "genes": genes}


def write_outputs(results, summary, out_dir):
    with open(os.path.join(out_dir, "benchmark_results.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "transcript", "hgvs_c", "expert_class", "tool_class", "agree_5class", "agree_3tier",
                    "expert_codes", "tool_codes", "points", "panel", "status", "note"])
        for r in results:
            w.writerow([r["gene"], r["transcript"], r["hgvs_c"], r["expert_class"], r.get("tool_class", ""),
                        r.get("tool_class") == r["expert_class"] if r["status"] == "ok" else "",
                        tier(r.get("tool_class", "")) == tier(r["expert_class"]) if r["status"] == "ok" else "",
                        ",".join(r["expert_codes"]), ",".join(r.get("tool_codes", [])), r.get("points", ""),
                        r.get("panel", ""), r["status"], r.get("note", "")])
    with open(os.path.join(out_dir, "benchmark_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
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
    tiers = "".join(f"<tr><td>{t}</td><td>{d['n']}</td><td>{pct(d['agree'] / d['n'] if d['n'] else None)}</td></tr>"
                    for t, d in summary["tier_agreement"].items())
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>CLAVAGEN benchmark</title>
<style>body{{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:28px;color:#1d2433;max-width:1100px}}
table{{border-collapse:collapse;margin:8px 0 22px}}td,th{{border:1px solid #e5e7eb;padding:6px 10px;text-align:left}}th{{background:#f3f4f6}}
.kpi{{display:flex;gap:14px;flex-wrap:wrap}}.kpi div{{background:#f4f6f9;border-radius:10px;padding:12px 16px;min-width:150px}}
.kpi b{{display:block;font-size:24px;color:#13294b}}.muted{{color:#6b7280;font-size:13px}}</style></head><body>
<h1>Benchmark against the ClinGen Evidence Repository</h1>
<p class="muted">{escape(summary['software'])} | {escape(summary['date'])} | automated ACMG engine only (ClinVar-based reported class OFF,
the variant's own ClinVar record excluded) | patient sequences built from RefSeq transcripts</p>
<div class="kpi"><div><b>{summary['analysed']}</b>variants analysed</div><div><b>{pct(summary['three_tier_agreement'])}</b>3-tier agreement
(P/LP, VUS, B/LB)</div><div><b>{pct(summary['exact_agreement'])}</b>exact 5-class agreement</div>
<div><b>{pct(summary['critical_discordance'])}</b>critical discordance (P/LP vs B/LB; n={summary['critical_n']})</div></div>
<p class="muted">Attempted {summary['attempted']}; not recovered {summary['not_recovered']}; errors {summary['errors']}.</p>
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
