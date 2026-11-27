"""Statistics for the ClinGen benchmark: confidence intervals, agreement coefficients, breakdowns.

Everything here is pure Python (no numpy/scipy), so the benchmark runs on a bare install.

* Proportions get Wilson score 95% intervals (Wilson 1927), which behave well for small n and
  for proportions near 0 or 1.
* Agreement beyond chance: Cohen's kappa on the 3-tier scale and linearly weighted kappa on the
  ordinal 5-class scale, with percentile bootstrap 95% intervals (fixed seed, reproducible).
"""

from __future__ import annotations

import math
import random

CLASSES = ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]
TIERS = ["P/LP", "VUS", "B/LB"]


def tier(c):
    return {"Pathogenic": "P/LP", "Likely pathogenic": "P/LP", "Uncertain significance": "VUS",
            "Likely benign": "B/LB", "Benign": "B/LB"}.get(c, "")


def wilson(k, n, z=1.959964):
    """Wilson score interval for k successes in n trials -> (p, low, high); None values if n == 0."""
    if not n:
        return None, None, None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return round(p, 4), round(max(0.0, c - h), 4), round(min(1.0, c + h), 4)


def prop(k, n):
    p, lo, hi = wilson(k, n)
    return {"k": k, "n": n, "p": p, "ci95": [lo, hi]}


def kappa(pairs, cats, weighted=False):
    """Cohen's kappa for (rater1, rater2) pairs; linear weights for ordinal categories if weighted."""
    n = len(pairs)
    if not n:
        return None
    idx = {c: i for i, c in enumerate(cats)}
    m = len(cats)
    obs = [[0] * m for _ in range(m)]
    for a, b in pairs:
        if a in idx and b in idx:
            obs[idx[a]][idx[b]] += 1
    tot = sum(map(sum, obs))
    if not tot:
        return None
    row = [sum(obs[i]) / tot for i in range(m)]
    col = [sum(obs[i][j] for i in range(m)) / tot for j in range(m)]

    def w(i, j):
        return 1 - abs(i - j) / (m - 1) if weighted else float(i == j)
    po = sum(w(i, j) * obs[i][j] / tot for i in range(m) for j in range(m))
    pe = sum(w(i, j) * row[i] * col[j] for i in range(m) for j in range(m))
    return None if pe >= 1 else (po - pe) / (1 - pe)


def kappa_ci(pairs, cats, weighted=False, n_boot=1000, seed=7):
    """Kappa with a percentile bootstrap 95% interval."""
    k = kappa(pairs, cats, weighted)
    if k is None or len(pairs) < 5:
        return {"kappa": None if k is None else round(k, 4), "ci95": [None, None], "n": len(pairs)}
    rnd = random.Random(seed)
    boots = []
    for _ in range(n_boot):
        s = [pairs[rnd.randrange(len(pairs))] for _ in pairs]
        b = kappa(s, cats, weighted)
        if b is not None:
            boots.append(b)
    boots.sort()
    lo = boots[int(0.025 * (len(boots) - 1))] if boots else None
    hi = boots[int(0.975 * (len(boots) - 1))] if boots else None
    return {"kappa": round(k, 4), "ci95": [None if lo is None else round(lo, 4), None if hi is None else round(hi, 4)],
            "n": len(pairs), "bootstrap": n_boot, "seed": seed}


def landis_koch(k):
    if k is None:
        return ""
    for t, lab in ((0.0, "poor"), (0.2, "slight"), (0.4, "fair"), (0.6, "moderate"), (0.8, "substantial")):
        if k <= t:
            return lab
    return "almost perfect"


def breakdown(ok, key, min_n=1):
    """Agreement by a grouping key (panel, gene, variant type, expert tier) with Wilson intervals."""
    groups = {}
    for r in ok:
        groups.setdefault(key(r) or "unknown", []).append(r)
    out = {}
    for g, rs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        if len(rs) < min_n:
            continue
        a3 = sum(1 for r in rs if tier(r["tool_class"]) == tier(r["expert_class"]))
        a5 = sum(1 for r in rs if r["tool_class"] == r["expert_class"])
        crit = sum(1 for r in rs if {tier(r["tool_class"]), tier(r["expert_class"])} == {"P/LP", "B/LB"})
        out[g] = {"n": len(rs), "three_tier": prop(a3, len(rs)), "exact": prop(a5, len(rs)),
                  "critical": prop(crit, len(rs)),
                  "expert_tiers": {t: sum(1 for r in rs if tier(r["expert_class"]) == t) for t in TIERS}}
    return out


def detection(ok, t):
    """Tool calls tier t vs the experts calling t: sensitivity and positive predictive value."""
    tp = sum(1 for r in ok if tier(r["tool_class"]) == t and tier(r["expert_class"]) == t)
    exp = sum(1 for r in ok if tier(r["expert_class"]) == t)
    called = sum(1 for r in ok if tier(r["tool_class"]) == t)
    return {"sensitivity": prop(tp, exp), "ppv": prop(tp, called)}


def direction(ok):
    """Among 5-class disagreements: is the tool more pathogenic (overcall) or more benign (undercall)?"""
    i = {c: k for k, c in enumerate(CLASSES)}
    over = sum(1 for r in ok if r["tool_class"] != r["expert_class"] and i[r["tool_class"]] < i[r["expert_class"]])
    under = sum(1 for r in ok if r["tool_class"] != r["expert_class"] and i[r["tool_class"]] > i[r["expert_class"]])
    t_over = sum(1 for r in ok if TIERS.index(tier(r["tool_class"])) < TIERS.index(tier(r["expert_class"])))
    t_under = sum(1 for r in ok if TIERS.index(tier(r["tool_class"])) > TIERS.index(tier(r["expert_class"])))
    return {"five_class": {"more_pathogenic": over, "more_benign": under},
            "three_tier": {"more_pathogenic": t_over, "more_benign": t_under}}


def coverage(ok):
    """How much automated evidence was available: population AF, REVEL for missense, transcript remapping."""
    mis = [r for r in ok if r.get("variant_type") == "missense"]
    return {"population_lookup": prop(sum(1 for r in ok if r.get("vep_resolved")), len(ok)),
            "population_af": prop(sum(1 for r in ok if r.get("population_af") is not None), len(ok)),
            "revel_missense": prop(sum(1 for r in mis if r.get("revel") is not None), len(mis)),
            "old_transcript_remapped": prop(sum(1 for r in ok if r.get("lookup_transcript")
                                                and r["lookup_transcript"] != r["transcript"]), len(ok))}


def noncomputable_use(ok):
    """How often the experts used criteria a sequence-based tool cannot evaluate."""
    from collections import Counter
    c = Counter(code for r in ok for code in r.get("expert_noncomputable", []))
    with_any = sum(1 for r in ok if r.get("expert_noncomputable"))
    return {"variants_with_any": prop(with_any, len(ok)), "by_code": dict(c.most_common())}


def supplemented(ok, seed=7):
    """Simulation: CLAVAGEN's automated criteria plus the clinical evidence the experts used (the
    criteria a sequence-based tool cannot compute: PS2-PS4, PM3, PM6, PP1, PP4, BS2-BS4, BP2, BP5...),
    as if a user had entered it with the guided clinical-evidence form. Upper bound of what entering
    clinical evidence can add; it does not test how users would describe that evidence."""
    rows = [r for r in ok if r.get("supplemented_class")]
    if not rows:
        return {}
    pairs3 = [(tier(r["expert_class"]), tier(r["supplemented_class"])) for r in rows]
    n = len(rows)
    with_ev = [r for r in rows if r.get("expert_noncomputable")]
    return {"n": n,
            "three_tier": prop(sum(1 for a, b in pairs3 if a == b), n),
            "exact": prop(sum(1 for r in rows if r["supplemented_class"] == r["expert_class"]), n),
            "critical_discordance": prop(sum(1 for a, b in pairs3 if {a, b} == {"P/LP", "B/LB"}), n),
            "kappa_three_tier": kappa_ci(pairs3, TIERS, seed=seed),
            "variants_with_clinical_evidence": len(with_ev),
            "three_tier_where_evidence_added": prop(sum(1 for r in with_ev if tier(r["supplemented_class"]) ==
                                                        tier(r["expert_class"])), len(with_ev)),
            "by_expert_tier": {t: prop(sum(1 for r in rows if tier(r["expert_class"]) == t and
                                           tier(r["supplemented_class"]) == t),
                                       sum(1 for r in rows if tier(r["expert_class"]) == t)) for t in TIERS}}


def publication_stats(ok, seed=7):
    pairs3 = [(tier(r["expert_class"]), tier(r["tool_class"])) for r in ok]
    pairs5 = [(r["expert_class"], r["tool_class"]) for r in ok]
    n = len(ok)
    return {
        "three_tier": prop(sum(1 for a, b in pairs3 if a == b), n),
        "exact": prop(sum(1 for a, b in pairs5 if a == b), n),
        "critical_discordance": prop(sum(1 for a, b in pairs3 if {a, b} == {"P/LP", "B/LB"}), n),
        "computable_three_tier": prop(sum(1 for r in ok if tier(r["tool_class"]) == tier(r["expert_computable_class"])), n),
        "computable_exact": prop(sum(1 for r in ok if r["tool_class"] == r["expert_computable_class"]), n),
        "kappa_three_tier": kappa_ci(pairs3, TIERS, seed=seed),
        "weighted_kappa_five_class": kappa_ci(pairs5, CLASSES, weighted=True, seed=seed),
        "detection": {t: detection(ok, t) for t in ("P/LP", "B/LB")},
        "direction": direction(ok),
        "by_expert_tier": breakdown(ok, lambda r: tier(r["expert_class"])),
        "by_panel": breakdown(ok, lambda r: r.get("panel")),
        "by_gene": breakdown(ok, lambda r: r.get("gene")),
        "by_variant_type": breakdown(ok, lambda r: r.get("variant_type")),
        "evidence_coverage": coverage(ok),
        "noncomputable_criteria": noncomputable_use(ok),
        "with_clinical_evidence": supplemented(ok, seed=seed),
    }
