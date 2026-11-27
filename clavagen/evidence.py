"""Guided clinical evidence: what the user observed -> ACMG/AMP criterion and strength.

Sequence data cannot show segregation, phase, functional assays, case counts or phenotype. Users who
have this information describe it here in plain terms (number of informative meioses, how often the
variant was seen in trans with a pathogenic variant, an assay's OddsPath...). CLAVAGEN converts it to
the criterion and strength with the published ClinGen rules, and shows the calculation in the report.

JSON (per variant key - c., p. or ``*`` for all variants of the patient):

{
  "c.20A>T": {
    "segregation":   {"informative_meioses": 5, "families": 2}            # or {"lod": 1.5}
    "non_segregation": true,                                              # BS4
    "in_trans": [{"partner": "Pathogenic", "phase": "confirmed"},        # PM3 (recessive)
                 {"partner": "Likely pathogenic", "phase": "unknown"}],
    "homozygous_probands": 1,                                             # PM3
    "benign_phase": true,                                                 # BP2
    "functional": {"result": "damaging", "oddspath": 20},                 # PS3 / BS3 (or "strength")
    "de_novo": [{"confirmed": true, "phenotype": "specific"}],            # PS2 / PM6
    "case_control": {"odds_ratio": 6.2, "ci_lower": 2.1},                 # PS4
    "affected_probands": {"count": 4, "strength": "Moderate"},            # PS4 (panel-specific counts)
    "phenotype_specific": true,                                           # PP4
    "healthy_adults": {"count": 2, "zygosity": "homozygous"},             # BS2
    "alternative_cause": true                                             # BP5
  }
}

References: Jarvik & Browning 2016 and Kelly et al. 2018 (segregation); ClinGen SVI PM3 v1.0 (2019) points;
Brnich et al. 2019 (functional evidence, OddsPath); ClinGen SVI de novo (PS2/PM6) points;
Richards et al. 2015 (PS4, BS2, BS4, BP2, BP5, PP4).
"""

from __future__ import annotations

import math

STRENGTHS = ("Supporting", "Moderate", "Strong", "VeryStrong")


def _by_points(points, table):
    """table: [(min_points, strength), ...] highest first."""
    for t, st in table:
        if points >= t - 1e-9:
            return st
    return None


DEFAULT_PP1_MEIOSES = {"Supporting": 3, "Moderate": 5, "Strong": 7}


def segregation(d, rules=None):
    """PP1 from informative meioses. Default: >=3 Supporting, >=5 Moderate, >=7 Strong (Kelly et al. 2018 /
    Jarvik & Browning 2016, as used by several ClinGen panels, e.g. PTEN, cardiomyopathy, glaucoma). Panels
    differ (e.g. thresholds that depend on the number of families), so a gene rule profile may set
    ``pp1_meioses``. A LOD score can be given instead: 3 meioses = LOD 0.9 (3 x log10 2)."""
    thr = dict(DEFAULT_PP1_MEIOSES, **((rules or {}).get("pp1_meioses") or {}))
    src = (rules or {}).get("source") if (rules or {}).get("pp1_meioses") else "default 3/5/7 meioses"
    n = d.get("informative_meioses")
    if n is None and d.get("lod") is not None:
        # LOD of k meioses = k x log10(2), conventionally rounded (3 -> 0.9, 5 -> 1.5, 7 -> 2.1)
        lod = float(d["lod"])
        n = max([0] + [k for k in range(1, 60) if lod >= round(k * math.log10(2), 1) - 1e-9])
        how = f"LOD {lod:.2f} (equivalent to {n} meioses)"
    elif n is not None:
        n = int(n)
        how = f"{n} informative meioses" + (f" in {d['families']} families" if d.get("families") else "")
    else:
        return []
    st = None
    for s_ in ("Strong", "Moderate", "Supporting"):
        if n >= thr[s_] - 1e-9:
            st = s_
            break
    rule = f"Supporting >= {thr['Supporting']}, Moderate >= {thr['Moderate']}, Strong >= {thr['Strong']} ({src})"
    if not st:
        return [("PP1", None, f"Co-segregation: {how} - fewer than needed for PP1 ({rule}).")]
    return [("PP1", st, f"Co-segregation with disease: {how}; {rule}.")]


def pm3(d):
    """PM3 points (ClinGen SVI recommendation v1.0): in trans with P confirmed 1.0, P phase unknown 0.5,
    LP confirmed 0.5, LP phase unknown 0.25; homozygous occurrence 0.5 each (max 1.0 in total).
    0.5 Supporting, 1.0 Moderate, 2.0 Strong, 4.0 VeryStrong."""
    pts, parts = 0.0, []
    for t in d.get("in_trans") or []:
        p = (t.get("partner") or "").lower()
        conf = (t.get("phase") or "unknown").lower() == "confirmed"
        if p.startswith("path"):
            x = 1.0 if conf else 0.5
        elif p.startswith("likely path") or p == "lp":
            x = 0.5 if conf else 0.25
        else:
            continue
        pts += x
        parts.append(f"{t.get('partner')} ({'in trans, confirmed' if conf else 'phase unknown'}) {x}")
    hom = int(d.get("homozygous_probands") or 0)
    if hom:
        x = min(1.0, 0.5 * hom)
        pts += x
        parts.append(f"{hom} homozygous proband(s) {x}")
    if not parts:
        return []
    st = _by_points(pts, [(4.0, "VeryStrong"), (2.0, "Strong"), (1.0, "Moderate"), (0.5, "Supporting")])
    why = (f"Recessive disorder, {pts:g} PM3 points: " + "; ".join(parts) +
           " (ClinGen SVI PM3: 0.5 Supporting, 1 Moderate, 2 Strong, 4 VeryStrong).")
    return [("PM3", st, why)] if st else [("PM3", None, why + " Below 0.5 point.")]


def functional(d):
    """PS3/BS3 (Brnich et al. 2019): OddsPath >18.7 Strong, >4.3 Moderate, >2.1 Supporting (PS3);
    <0.053 Strong, <0.23 Moderate, <0.48 Supporting (BS3). Without an OddsPath the user's strength is used
    (default Strong for a well-established assay, per ACMG 2015)."""
    res = (d.get("result") or "").lower()
    if res not in ("damaging", "abnormal", "normal", "benign", "no effect"):
        return []
    path = res in ("damaging", "abnormal")
    code = "PS3" if path else "BS3"
    op = d.get("oddspath")
    if op is not None:
        op = float(op)
        st = (_by_points(op, [(18.7, "Strong"), (4.3, "Moderate"), (2.1, "Supporting")]) if path else
              next((s for t, s in ((0.053, "Strong"), (0.23, "Moderate"), (0.48, "Supporting")) if op < t), None))
        why = f"Functional assay ({res}), OddsPath {op:g} (Brnich et al. 2019 thresholds)."
        return [(code, st, why if st else why + " Assay not sufficiently validated for this criterion.")]
    st = d.get("strength") or "Strong"
    return [(code, st, f"Well-established functional assay shows a {res} effect"
                       f"{' (' + d['assay'] + ')' if d.get('assay') else ''}; strength {st} as entered.")]


def de_novo(probands):
    """PS2/PM6 points (ClinGen SVI de novo): phenotype highly specific - confirmed 2, assumed 1;
    consistent but not specific - confirmed 1, assumed 0.5; ... Points 0.5 Supporting, 1 Moderate,
    2 Strong, 4 VeryStrong. Confirmed de novo (maternity and paternity) -> PS2, otherwise PM6."""
    if isinstance(probands, dict):
        probands = [probands]
    pts, low, conf_any, parts = 0.0, 0.0, False, []
    for p in probands or []:
        conf = bool(p.get("confirmed"))
        ph = (p.get("phenotype") or "consistent").lower()
        base = {"specific": 2.0, "consistent": 1.0}.get(ph, 0.5 if ph.startswith("hetero") else 1.0)
        x = base if conf else base / 2
        if base == 0.5:      # low-specificity phenotypes contribute at most 1 point in total
            x = max(0.0, min(x, 1.0 - low))
            low += x
        pts += x
        conf_any |= conf
        parts.append(f"{'confirmed' if conf else 'assumed'} de novo, phenotype {ph}: {x:g}")
    if not parts:
        return []
    st = _by_points(pts, [(4.0, "VeryStrong"), (2.0, "Strong"), (1.0, "Moderate"), (0.5, "Supporting")])
    code = "PS2" if conf_any else "PM6"
    return [(code, st, f"De novo occurrence, {pts:g} points: " + "; ".join(parts) +
             " (ClinGen SVI de novo points: 0.5 Supporting, 1 Moderate, 2 Strong, 4 VeryStrong).")]


def build(d, rules=None):
    """Structured evidence for one variant -> [(code, strength or None, rationale)].

    A strength of None means the evidence was considered but does not reach the criterion; it is reported
    as 'not met' with the calculation."""
    out = []
    if d.get("segregation"):
        out += segregation(d["segregation"], rules)
    if d.get("non_segregation"):
        out.append(("BS4", "Strong", "Lack of segregation in affected family members (entered by user)."))
    if d.get("in_trans") or d.get("homozygous_probands"):
        out += pm3(d)
    if d.get("benign_phase"):
        out.append(("BP2", "Supporting", "Observed in trans with a pathogenic variant in a dominant disorder, or in cis "
                                         "with a pathogenic variant (entered by user)."))
    if d.get("functional"):
        out += functional(d["functional"])
    if d.get("de_novo"):
        out += de_novo(d["de_novo"])
    cc = d.get("case_control") or {}
    if cc.get("odds_ratio") is not None:
        ok = float(cc["odds_ratio"]) >= 5 and float(cc.get("ci_lower", 0)) > 1
        out.append(("PS4", "Strong" if ok else None,
                    f"Case-control odds ratio {cc['odds_ratio']} (95% CI lower {cc.get('ci_lower', '?')}); ACMG PS4 "
                    "needs OR >= 5 with a CI not including 1."))
    ap = d.get("affected_probands") or {}
    if ap.get("count"):
        st = ap.get("strength") or "Supporting"
        out.append(("PS4", st, f"{ap['count']} unrelated affected probands (strength {st} per the gene's panel "
                               "specification, as entered)."))
    if d.get("phenotype_specific"):
        out.append(("PP4", "Supporting", "Patient phenotype or family history highly specific for a disease with a "
                                         "single genetic aetiology (entered by user)."))
    ha = d.get("healthy_adults") or {}
    if ha.get("count"):
        st = ha.get("strength") or "Strong"
        out.append(("BS2", st, f"Observed in {ha['count']} healthy adult(s) ({ha.get('zygosity', 'genotype as entered')}) "
                               "for a fully penetrant, early-onset disorder."))
    if d.get("alternative_cause"):
        out.append(("BP5", "Supporting", "An alternative molecular cause of the disease was found (entered by user)."))
    # one entry per code: the strongest wins
    rank = {None: -1, "Supporting": 0, "Moderate": 1, "Strong": 2, "VeryStrong": 3, "StandAlone": 4}
    best = {}
    for code, st, why in out:
        if code not in best or rank[st] > rank[best[code][0]]:
            best[code] = (st, why)
    return [(c, s, w) for c, (s, w) in best.items()]


def from_json(obj, rules=None):
    """{key: evidence dict} -> {key: [(code, strength, rationale)]}"""
    return {k: build(v, rules) for k, v in (obj or {}).items() if isinstance(v, dict)}
