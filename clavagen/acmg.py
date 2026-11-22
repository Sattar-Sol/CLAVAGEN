"""ACMG/AMP 2015 variant classification engine with ClinGen SVI refinements.

References
----------
* Richards S. et al. Genet Med 2015;17:405-424 (ACMG/AMP standards; combining rules, Table 5).
* Abou Tayoun A.N. et al. Hum Mutat 2018;39:1517 (PVS1 decision tree).
* ClinGen SVI recommendation for PM2 (2020): apply at Supporting strength.
* Tavtigian S.V. et al. Hum Mutat 2018; 2020 (Bayesian / points framework).
* Pejaver V. et al. AJHG 2022;109:2163 (calibrated REVEL thresholds for PP3/BP4).
* Walker L.C. et al. AJHG 2023 (SpliceAI / splicing evidence; BP7 for deep intronic).
PP5 and BP6 are not used, following ClinGen SVI (2018).

Only criteria that can be assessed from sequence + annotation are automated;
all others are accepted as user-supplied evidence. Heuristic (non-calibrated)
evidence is flagged and can be disabled with ``strict=True``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .annotation import (insilico_lookup, is_pathogenic_class, normalize_hgvs, parse_missense,
                         population_lookup)
from .consequence import LOF_CONSEQUENCES
from .genetic_code import AA_CLASS

STRENGTHS = ["StandAlone", "VeryStrong", "Strong", "Moderate", "Supporting"]
POINTS = {"VeryStrong": 8, "Strong": 4, "Moderate": 2, "Supporting": 1, "StandAlone": 8}

CRITERIA = {
    # code: (direction, default strength, short description)
    "PVS1": ("P", "VeryStrong", "Null variant in a gene where loss of function is a known disease mechanism"),
    "PS1": ("P", "Strong", "Same amino-acid change as an established pathogenic variant"),
    "PS2": ("P", "Strong", "De novo (maternity and paternity confirmed)"),
    "PS3": ("P", "Strong", "Well-established functional studies show a damaging effect"),
    "PS4": ("P", "Strong", "Prevalence in affected individuals significantly increased vs controls"),
    "PM1": ("P", "Moderate", "Located in a mutational hot spot / critical functional domain"),
    "PM2": ("P", "Moderate", "Absent / extremely rare in population databases"),
    "PM3": ("P", "Moderate", "For recessive disorders, detected in trans with a pathogenic variant"),
    "PM4": ("P", "Moderate", "Protein length change (in-frame indel in non-repeat region or stop-loss)"),
    "PM5": ("P", "Moderate", "Novel missense at a residue where a different pathogenic missense is known"),
    "PM6": ("P", "Moderate", "Assumed de novo (without confirmation of parentage)"),
    "PP1": ("P", "Supporting", "Co-segregation with disease in multiple affected family members"),
    "PP2": ("P", "Supporting", "Missense in a gene with low benign missense rate where missense is a mechanism"),
    "PP3": ("P", "Supporting", "Computational evidence supports a deleterious effect"),
    "PP4": ("P", "Supporting", "Patient phenotype / family history highly specific for the gene"),
    "BA1": ("B", "StandAlone", "Allele frequency > 5% (or gene-specific BA1 threshold)"),
    "BS1": ("B", "Strong", "Allele frequency greater than expected for the disorder"),
    "BS2": ("B", "Strong", "Observed in healthy adult(s) for a fully penetrant early-onset disorder"),
    "BS3": ("B", "Strong", "Well-established functional studies show no damaging effect"),
    "BS4": ("B", "Strong", "Lack of segregation in affected family members"),
    "BP1": ("B", "Supporting", "Missense in a gene where only truncating variants cause disease"),
    "BP2": ("B", "Supporting", "In trans with a pathogenic variant (dominant gene) or in cis with a pathogenic variant"),
    "BP3": ("B", "Supporting", "In-frame indel in a repetitive region without known function"),
    "BP4": ("B", "Supporting", "Computational evidence suggests no impact"),
    "BP5": ("B", "Supporting", "Found in a case with an alternate molecular basis for disease"),
    "BP7": ("B", "Supporting", "Synonymous/deep-intronic with no predicted splice impact"),
}

# Pejaver et al. 2022 REVEL calibration
REVEL_PP3 = [(0.932, "Strong"), (0.773, "Moderate"), (0.644, "Supporting")]
REVEL_BP4 = [(0.003, "VeryStrong"), (0.016, "Strong"), (0.183, "Moderate"), (0.290, "Supporting")]

DEFAULT_THRESHOLDS = {
    "ba1": 0.05,
    "bs1": {"AD": 0.001, "AR": 0.01, "XL": 0.001, "default": 0.005},
    "pm2": {"AD": 0.00001, "AR": 0.0005, "XL": 0.00001, "default": 0.0001},
}


@dataclass
class Criterion:
    code: str
    status: str            # met | not_met | not_evaluated
    strength: str = ""
    rationale: str = ""
    source: str = "automated"  # automated | heuristic | user | database

    @property
    def direction(self):
        return CRITERIA.get(self.code, ("P", "", ""))[0]

    @property
    def label(self):
        default = CRITERIA.get(self.code, ("", "", ""))[1]
        if not self.strength or self.strength == default:
            return self.code
        return f"{self.code}_{self.strength}"

    @property
    def points(self):
        if self.status != "met":
            return 0
        p = POINTS.get(self.strength, 0)
        return p if self.direction == "P" else -p

    def to_dict(self):
        d = asdict(self)
        d.update(label=self.label, direction=self.direction, points=self.points,
                 description=CRITERIA.get(self.code, ("", "", ""))[2])
        return d


def _inh_key(inh):
    inh = (inh or "").upper()
    if "AR" in inh or "RECESS" in inh:
        return "AR"
    if inh.startswith("X") or "X-" in inh:
        return "XL"
    if "AD" in inh or "DOMIN" in inh:
        return "AD"
    return "default"


def _threshold(ctx, name):
    th = (ctx.get("thresholds") or {}).get(name)
    if th is not None:
        return float(th)
    d = DEFAULT_THRESHOLDS[name]
    if isinstance(d, dict):
        return d.get(_inh_key(ctx.get("inheritance")), d["default"])
    return d


def parse_user_evidence(items):
    """Parse ['c.20A>T=PS3,PP1_Strong', '*=PP4', 'p.Glu7Val=-PM2'] -> {key: [(code, strength|None, remove)]}."""
    ev = {}
    for it in items or []:
        if "=" not in it:
            continue
        key, codes = it.split("=", 1)
        key = "*" if key.strip() == "*" else normalize_hgvs(key.strip())
        for c in codes.replace(";", ",").split(","):
            c = c.strip()
            if not c:
                continue
            remove = c.startswith("-")
            c = c.lstrip("-+")
            code, _, strength = c.partition("_")
            code = code.upper()
            strength = {"VERYSTRONG": "VeryStrong", "VSTRONG": "VeryStrong", "STRONG": "Strong",
                        "MODERATE": "Moderate", "SUPPORTING": "Supporting", "STANDALONE": "StandAlone"}.get(
                strength.upper(), None) if strength else None
            if code in CRITERIA:
                ev.setdefault(key, []).append((code, strength, remove))
    return ev


# ------------------------------------------------------------------ PVS1
def _pvs1(v, ctx):
    cons = v["consequence"]
    plen = v.get("protein_length") or 0
    junc = v.get("junctions_c") or []
    if cons == "start_lost":
        return "Moderate", ("Start-loss: PVS1_Moderate per ClinGen PVS1 tree (check for an in-frame "
                            "downstream Met and pathogenic variants upstream of it).")
    if cons in ("splice_donor_variant", "splice_acceptor_variant", "exon_loss_variant"):
        return "Strong", ("Canonical +/-1,2 splice site or exon-boundary change. Default PVS1_Strong; may be "
                          "upgraded to VeryStrong if exon skipping is out-of-frame and NMD is predicted, or "
                          "downgraded if the exon is in-frame and not critical. Confirm with RNA studies/SpliceAI.")
    ptc = v.get("ptc_pos")
    if ptc is None and v.get("aa_pos"):
        ptc = v["aa_pos"]
    ptc_c = (ptc - 1) * 3 + 1 if ptc else None
    if junc and ptc_c:
        last_j = max(junc)
        if ptc_c <= last_j - 50:
            return "VeryStrong", (f"Premature stop at codon {ptc} (c.{ptc_c}) lies >50 nt upstream of the last "
                                  f"exon-exon junction (c.{last_j}/c.{last_j + 1}); NMD predicted.")
        frac = (plen - ptc + 1) / plen if plen else 0
        msg = (f"PTC at codon {ptc} is in the last exon or within 50 nt of the last junction: NMD escape "
               f"predicted; ~{frac * 100:.0f}% of the protein removed")
        crit_dom = [d for d in v.get("domains_downstream", []) if d.get("critical")]
        if frac > 0.10 or crit_dom:
            return "Strong", msg + (" (critical region removed)." if crit_dom else " (>10%).")
        return "Moderate", msg + " (<=10%)."
    frac = (plen - ptc + 1) / plen if (plen and ptc) else 1.0
    msg = ("Exon structure unknown - NMD cannot be predicted (supply 'exons' or 'exon_boundaries_c'). "
           f"Truncation removes ~{frac * 100:.0f}% of the protein.")
    if frac > 0.10:
        return "Strong", msg
    return "Moderate", msg


# ------------------------------------------------------------------ main
def evaluate(v, ctx, user_ev=None, strict=False, patient_phenotypes=None, family_ev=None):
    """Evaluate ACMG/AMP criteria for one annotated variant (dict). Returns list[Criterion]."""
    C = []
    cons = v["consequence"]
    keys = [k for k in (v.get("hgvs_c"), v.get("hgvs_p"), v.get("hgvs_ref")) if k]
    nkeys = {normalize_hgvs(k) for k in keys}
    is_missense = cons == "missense_variant"
    inframe = cons in ("inframe_deletion", "inframe_insertion", "protein_altering_variant")

    # ---------------- PVS1
    if cons in LOF_CONSEQUENCES:
        lof = ctx.get("lof_mechanism")
        if lof is True:
            st, why = _pvs1(v, ctx)
            C.append(Criterion("PVS1", "met", st, why))
        elif lof is False:
            C.append(Criterion("PVS1", "not_met", "", "Null variant, but LOF is not an established mechanism for this gene."))
        else:
            C.append(Criterion("PVS1", "not_evaluated", "",
                               "Null variant; gene LOF mechanism unknown - set 'lof_mechanism' in annotation."))

    # ---------------- known variants: PS1 / PM5 / database assertions
    known = ctx.get("known_variants") or []
    v_mis = (v.get("aa_pos"), v.get("ref_aa"), v.get("alt_aa")) if is_missense else None
    nuc_known = bool(v.get("hgvs_c"))                      # DNA input, or uniquely inferred from protein
    cands = {normalize_hgvs(c) for c in v.get("hgvs_c_candidates", [])}
    pkey = normalize_hgvs(v.get("hgvs_p", ""))
    ps1 = pm5 = None
    matches = []
    for kv in known:
        kc, kp = normalize_hgvs(kv.get("hgvs_c", "")), normalize_hgvs(kv.get("hgvs_p", ""))
        exact = bool(kc and (kc in nkeys or kc in cands))
        prot_only = bool(not exact and kp and kp == pkey and (not kc or not nuc_known))
        if exact or prot_only:
            matches.append((exact, kv))
            star = f", {kv['stars']}\u2605 {kv.get('review_status', '')}" if kv.get("stars") is not None and kv.get("review_status") else ""
            tag = "" if exact else " [same protein change; nucleotide not confirmed]"
            v.setdefault("database_assertions", []).append(
                f"{kv.get('classification', '?')} ({kv.get('source', 'user database')}{star}){tag}")
        if not is_missense or not is_pathogenic_class(kv.get("classification")):
            continue
        km = parse_missense(kv.get("hgvs_p", ""))
        if not km:
            continue
        if km == v_mis and not exact:
            # same amino-acid change, different (or unconfirmed) nucleotide change
            if not ps1 or (kv.get("stars", 0) or 0) > (ps1.get("stars", 0) or 0):
                ps1 = kv
        elif km[0] == v_mis[0] and km[1] == v_mis[1] and km[2] != v_mis[2] and km[2] != "*":
            if not pm5 or (kv.get("stars", 0) or 0) > (pm5.get("stars", 0) or 0):
                pm5 = kv
    # best database record for the variant itself (exact match preferred, then review stars)
    if matches:
        matches.sort(key=lambda m: (m[0], not str(m[1].get("source", "")).startswith("ClinVar"),
                                    m[1].get("stars", 0) or 0), reverse=True)
        best = dict(matches[0][1])
        best["match"] = "exact" if matches[0][0] else "protein"
        v["clinvar_match"] = best
    if is_missense:
        if ps1:
            if nuc_known:
                why = (f"Same amino-acid change as {ps1.get('hgvs_c', '')} {ps1.get('hgvs_p', '')} "
                       f"({ps1.get('classification')}, {ps1.get('source', 'database')}) with a different nucleotide "
                       "change; confirm no splicing effect of either change.")
            else:
                why = (f"Same amino-acid change as established {ps1.get('classification')} variant "
                       f"{ps1.get('hgvs_c', '')} {ps1.get('hgvs_p', '')} ({ps1.get('source', 'database')}). Protein-level "
                       "input: the nucleotide change is not confirmed - if it is the identical variant, use the "
                       "database classification instead of PS1.")
            C.append(Criterion("PS1", "met", "Strong", why, "database"))
        if pm5 and not ps1:
            C.append(Criterion("PM5", "met", "Moderate",
                               f"Different pathogenic missense at the same residue: {pm5.get('hgvs_p')} "
                               f"({pm5.get('classification')}, {pm5.get('source', 'database')}).", "database"))
        elif pm5 and ps1:
            C.append(Criterion("PM5", "not_met", "", f"PS1 applied; PM5 ({pm5.get('hgvs_p')}) not combined with PS1."))
        if not known:
            C.append(Criterion("PS1", "not_evaluated", "", "No known-variant data (online ClinVar search off/unavailable "
                                                           "and no manual list)."))

    # ---------------- PM1 (domain / hotspot)
    if is_missense or inframe:
        hs = set(int(x) for x in ctx.get("hotspots", []) or [])
        crit = [d for d in v.get("domains", []) if d.get("critical")]
        sites = set(ctx.get("critical_sites") or [])
        if v.get("aa_pos") in hs:
            C.append(Criterion("PM1", "met", "Moderate",
                               f"Residue {v['aa_pos']} is a mutational hotspot (annotated, or >=2 distinct P/LP "
                               "missense changes at this residue in ClinVar)."))
        elif is_missense and v.get("aa_pos") in sites:
            C.append(Criterion("PM1", "met", "Moderate",
                               f"Residue {v['aa_pos']} is an active/binding/functional site (UniProt)."))
        elif crit:
            C.append(Criterion("PM1", "met", "Moderate",
                               f"Located in critical functional domain: {', '.join(d['name'] for d in crit)}."))
        elif v.get("domains"):
            C.append(Criterion("PM1", "not_met", "",
                               f"In domain {', '.join(d['name'] for d in v['domains'])} (not flagged critical/hotspot)."))

    # ---------------- population: BA1 / BS1 / PM2
    pop = population_lookup(ctx, keys)
    has_pop_data = bool(ctx.get("population")) or bool(ctx.get("population_dataset"))
    ba1_t, bs1_t, pm2_t = _threshold(ctx, "ba1"), _threshold(ctx, "bs1"), _threshold(ctx, "pm2")
    ds = ctx.get("population_dataset") or "provided population data"
    if pop is not None:
        af = float(pop.get("af", 0))
        src = pop.get("source", ds)
        exc = bool(nkeys & ctx.get("ba1_exceptions", set()))
        if af > ba1_t and not exc:
            C.append(Criterion("BA1", "met", "StandAlone", f"AF {af:.4g} > BA1 threshold {ba1_t} ({src})."))
        elif af > ba1_t and exc:
            C.append(Criterion("BA1", "not_met", "", f"AF {af:.4g} > {ba1_t} but variant is on the ClinGen BA1 exception list."))
        elif af > bs1_t:
            C.append(Criterion("BS1", "met", "Strong", f"AF {af:.4g} > disorder-specific BS1 threshold {bs1_t} ({src})."))
        elif af <= pm2_t:
            C.append(Criterion("PM2", "met", "Supporting",
                               f"AF {af:.4g} <= PM2 threshold {pm2_t} ({src}); PM2 applied at Supporting (ClinGen SVI 2020)."))
        else:
            C.append(Criterion("PM2", "not_met", "", f"AF {af:.4g} above PM2 threshold {pm2_t} but below BS1 {bs1_t}."))
    elif has_pop_data and v.get("mol_type") == "protein" and not v.get("hgvs_c"):
        C.append(Criterion("PM2", "not_evaluated", "",
                           "Protein-level input: population frequency requires the nucleotide change "
                           "(supply a p.-keyed entry in 'population' or analyse the DNA sequence)."))
    elif has_pop_data:
        C.append(Criterion("PM2", "met", "Supporting",
                           f"Absent from {ds}; PM2 applied at Supporting (ClinGen SVI 2020). "
                           "Confirm adequate coverage at this position."))
    else:
        C.append(Criterion("PM2", "not_evaluated", "",
                           "No population data supplied - query gnomAD and add to annotation 'population'."))

    # ---------------- PM4 / BP3
    rep = (v.get("repeat_copies") or 0) >= 3 or (v.get("repeat_copies", 0) >= 2 and len(v.get("repeat_unit", "")) >= 3)
    if inframe:
        if rep and not any(d.get("critical") for d in v.get("domains", [])):
            C.append(Criterion("BP3", "met", "Supporting",
                               f"In-frame change in a repetitive region ({v.get('repeat_unit')} x{v.get('repeat_copies')}) "
                               "outside annotated critical domains."))
        elif not rep:
            C.append(Criterion("PM4", "met", "Moderate", "In-frame indel/delins in a non-repetitive region alters protein length."))
    if cons == "stop_lost":
        C.append(Criterion("PM4", "met", "Moderate", "Stop-loss extends the protein."))

    # ---------------- PP2 / BP1
    if is_missense:
        if ctx.get("missense_constrained"):
            C.append(Criterion("PP2", "met", "Supporting", "Gene annotated as missense-constrained with missense as a disease mechanism."))
        if ctx.get("truncating_only"):
            C.append(Criterion("BP1", "met", "Supporting", "Gene annotated: truncating variants are the primary disease mechanism."))

    # ---------------- PP3 / BP4 (in silico)
    ins = insilico_lookup(ctx, keys) or {}
    revel = ins.get("REVEL") if isinstance(ins, dict) else None
    spliceai = ins.get("SpliceAI") if isinstance(ins, dict) else None
    pp3 = bp4 = None
    if is_missense and revel is not None:
        r = float(revel)
        for t, st in REVEL_PP3:
            if r >= t:
                pp3 = Criterion("PP3", "met", st, f"REVEL {r:.3f} >= {t} (Pejaver 2022 calibration).")
                break
        if not pp3:
            for t, st in REVEL_BP4:
                if r <= t:
                    bp4 = Criterion("BP4", "met", st, f"REVEL {r:.3f} <= {t} (Pejaver 2022 calibration).")
                    break
        if not pp3 and not bp4:
            C.append(Criterion("PP3", "not_met", "", f"REVEL {r:.3f} in indeterminate range (0.290-0.644)."))
    if spliceai is not None and not pp3:
        s = float(spliceai)
        if s >= 0.2:
            pp3 = Criterion("PP3", "met", "Supporting", f"SpliceAI max delta {s:.2f} >= 0.2 (predicted splice impact).")
        elif s <= 0.1 and not is_missense and not bp4:
            bp4 = Criterion("BP4", "met", "Supporting", f"SpliceAI max delta {s:.2f} <= 0.1 (no predicted splice impact).")
    if is_missense and revel is None and not strict:
        h = heuristic_missense(v)
        if h["call"] == "damaging":
            pp3 = Criterion("PP3", "met", "Supporting", "Heuristic (not calibrated): " + h["why"], "heuristic")
        elif h["call"] == "tolerated":
            bp4 = Criterion("BP4", "met", "Supporting", "Heuristic (not calibrated): " + h["why"], "heuristic")
        else:
            C.append(Criterion("PP3", "not_met", "", "Heuristic in-silico evidence indeterminate: " + h["why"], "heuristic"))
    elif is_missense and revel is None:
        C.append(Criterion("PP3", "not_evaluated", "", "No REVEL score supplied (strict mode: heuristics disabled)."))
    if pp3:
        C.append(pp3)
    if bp4:
        C.append(bp4)

    # ---------------- BP7 (synonymous / deep intronic)
    if cons == "synonymous_variant" or (cons == "intron_variant" and v.get("intron_offset") is not None
                                         and (v["intron_offset"] >= 7 or v["intron_offset"] <= -21)):
        if v.get("splice_region"):
            C.append(Criterion("BP7", "not_met", "", "Within splice region - splicing impact cannot be excluded."))
        elif spliceai is not None and float(spliceai) > 0.1:
            C.append(Criterion("BP7", "not_met", "", f"SpliceAI {spliceai} > 0.1."))
        else:
            src = "automated" if spliceai is not None else "heuristic"
            C.append(Criterion("BP7", "met", "Supporting",
                               ("Synonymous" if cons == "synonymous_variant" else "Deep intronic") +
                               " variant outside the splice region" +
                               (f"; SpliceAI {spliceai}" if spliceai is not None else
                                " (no splice predictor supplied - confirm with SpliceAI)") + ".", src))

    # ---------------- PP4 (phenotype)
    pp4_eligible = cons not in ("synonymous_variant", "stop_retained_variant", "intron_variant",
                                "5_prime_UTR_variant", "3_prime_UTR_variant", "upstream_gene_variant",
                                "downstream_gene_variant", "non_coding_variant")
    if patient_phenotypes and pp4_eligible:
        kw = [k.lower() for k in (ctx.get("phenotype_keywords") or [])]
        feats = [f.lower() for d in ctx.get("diseases", []) for f in d.get("features", [])]
        hits = sorted({p for p in patient_phenotypes if any(k in p.lower() for k in kw) or
                       any(p.lower() in f or f in p.lower() for f in feats)})
        if hits and not strict:
            C.append(Criterion("PP4", "met", "Supporting",
                               f"Patient phenotype matches the gene-disease spectrum ({'; '.join(hits)}). "
                               "Apply only if phenotype is highly specific for a single genetic aetiology.", "heuristic"))
        elif not hits:
            C.append(Criterion("PP4", "not_met", "", "Supplied phenotype does not match the annotated gene-disease spectrum."))

    # ---------------- PM3 hint
    if _inh_key(ctx.get("inheritance")) == "AR":
        C.append(Criterion("PM3", "not_evaluated", "",
                           "AR gene: evaluate PM3 if a pathogenic variant is confirmed in trans (phasing required)."))

    # ---------------- SVI cap: PP3 + PM1 combined <= Strong (4 points)
    met = {c.code: c for c in C if c.status == "met"}
    if "PP3" in met and "PM1" in met and POINTS[met["PP3"].strength] + POINTS[met["PM1"].strength] > 4:
        met["PM1"].strength = "Supporting"
        met["PM1"].rationale += " [Downgraded: combined PP3+PM1 capped at Strong (ClinGen SVI).]"

    # ---------------- family evidence (de novo, phase, segregation)
    for code, strength, why in (family_ev or []):
        C = [c for c in C if c.code != code]
        C.append(Criterion(code, "met", strength, why, "family"))

    # ---------------- user evidence (overrides automated)
    ue = user_ev or {}
    for key in list(nkeys) + ["*"]:
        for code, strength, remove in ue.get(key, []):
            C = [c for c in C if c.code != code]
            if remove:
                C.append(Criterion(code, "not_met", "", "Removed by reviewer.", "user"))
            else:
                C.append(Criterion(code, "met", strength or CRITERIA[code][1], "Supplied by reviewer/clinician.", "user"))
    return C


def heuristic_missense(v):
    """Simple, transparent, non-calibrated in-silico call based on BLOSUM62, Grantham,
    physico-chemical class change and (optionally) ortholog conservation."""
    pts, why = 0, []
    b, g = v.get("blosum62"), v.get("grantham")
    r, a = v.get("ref_aa", ""), v.get("alt_aa", "")
    if b is not None:
        if b <= -2:
            pts += 1
            why.append(f"BLOSUM62 {b}")
        elif b >= 1:
            pts -= 1
            why.append(f"BLOSUM62 +{b}")
    if g is not None:
        if g >= 100:
            pts += 1
            why.append(f"Grantham {g}")
        elif g < 60:
            pts -= 1
            why.append(f"Grantham {g}")
    rc, ac = AA_CLASS.get(r), AA_CLASS.get(a)
    if rc and ac and rc != ac:
        if {rc, ac} == {"positive", "negative"} or "special" in (rc, ac) or rc in ("positive", "negative"):
            pts += 1
        why.append(f"{rc}->{ac}")
    cons = v.get("conservation")
    if cons and cons.get("n_orthologs", 0) >= 3:
        if cons.get("ref_identity", 0) >= 0.99 and not cons.get("alt_observed"):
            pts += 1
            why.append(f"ref residue invariant in {cons['n_orthologs']} orthologs")
        elif cons.get("alt_observed"):
            pts -= 2
            why.append("alt residue observed in orthologs")
    need = 4 if (cons and cons.get("n_orthologs", 0) >= 3) else 3
    call = "damaging" if pts >= need else ("tolerated" if pts <= -2 else "indeterminate")
    return {"call": call, "points": pts, "why": ", ".join(why) or "no informative features"}


# ------------------------------------------------------------------ combine
def combine(criteria):
    """Apply ACMG/AMP combining rules (with SVI PVS1+PP->LP) and the points system."""
    met = [c for c in criteria if c.status == "met"]
    cnt = {"PVS": 0, "PS": 0, "PM": 0, "PP": 0, "BA": 0, "BS": 0, "BP": 0}
    for c in met:
        if c.direction == "P":
            k = {"VeryStrong": "PVS", "Strong": "PS", "Moderate": "PM", "Supporting": "PP"}.get(c.strength, "PP")
            cnt[k] += 1
        else:
            if c.strength == "StandAlone":
                cnt["BA"] += 1
            elif c.strength in ("VeryStrong", "Strong"):
                cnt["BS"] += 1
            elif c.strength == "Moderate":
                cnt["BP"] += 2
            else:
                cnt["BP"] += 1
    PVS, PS, PM, PP = cnt["PVS"], cnt["PS"], cnt["PM"], cnt["PP"]
    path, rule = None, ""
    if PVS >= 1 and (PS >= 1 or PM >= 2 or (PM == 1 and PP >= 1) or PP >= 2):
        path, rule = "Pathogenic", "P(i): PVS1 + (>=1 PS | >=2 PM | PM+PP | >=2 PP)"
    elif PS >= 2:
        path, rule = "Pathogenic", "P(ii): >=2 PS"
    elif PS == 1 and (PM >= 3 or (PM == 2 and PP >= 2) or (PM == 1 and PP >= 4)):
        path, rule = "Pathogenic", "P(iii): PS + (>=3 PM | 2 PM + >=2 PP | PM + >=4 PP)"
    elif PVS >= 1 and PM == 1:
        path, rule = "Likely pathogenic", "LP(i): PVS1 + PM"
    elif PVS >= 1 and PP >= 1:
        path, rule = "Likely pathogenic", "LP (ClinGen SVI): PVS1 + PP (e.g. PM2_Supporting)"
    elif PS == 1 and 1 <= PM <= 2:
        path, rule = "Likely pathogenic", "LP(ii): PS + 1-2 PM"
    elif PS == 1 and PP >= 2:
        path, rule = "Likely pathogenic", "LP(iii): PS + >=2 PP"
    elif PM >= 3:
        path, rule = "Likely pathogenic", "LP(iv): >=3 PM"
    elif PM == 2 and PP >= 2:
        path, rule = "Likely pathogenic", "LP(v): 2 PM + >=2 PP"
    elif PM == 1 and PP >= 4:
        path, rule = "Likely pathogenic", "LP(vi): PM + >=4 PP"
    ben, brule = None, ""
    if cnt["BA"] >= 1:
        ben, brule = "Benign", "B(i): BA1"
    elif cnt["BS"] >= 2:
        ben, brule = "Benign", "B(ii): >=2 BS"
    elif cnt["BS"] == 1 and cnt["BP"] >= 1:
        ben, brule = "Likely benign", "LB(i): BS + BP"
    elif cnt["BP"] >= 2:
        ben, brule = "Likely benign", "LB(ii): >=2 BP"

    if cnt["BA"] >= 1:
        rules_class, used = "Benign", brule
    elif path and ben:
        rules_class, used = "Uncertain significance", f"Conflicting evidence ({rule} vs {brule})"
    elif path:
        rules_class, used = path, rule
    elif ben:
        rules_class, used = ben, brule
    else:
        rules_class, used = "Uncertain significance", "Criteria for other classifications not met"

    points = sum(c.points for c in met if c.code != "BA1")
    if points >= 10:
        pclass = "Pathogenic"
    elif points >= 6:
        pclass = "Likely pathogenic"
    elif points >= 0:
        pclass = "Uncertain significance"
    elif points >= -6:
        pclass = "Likely benign"
    else:
        pclass = "Benign"
    if cnt["BA"]:
        pclass = "Benign"
    sub = ""
    if rules_class == "Uncertain significance":
        sub = "VUS-high" if points >= 4 else ("VUS-mid" if points >= 2 else "VUS-low")
    return {
        "classification": rules_class,
        "rule": used,
        "points": points,
        "points_classification": pclass,
        "vus_subtier": sub,
        "concordant": rules_class == pclass,
        "criteria_met": [c.label for c in met],
        "counts": cnt,
    }


CLASS_ABBR = {"Pathogenic": "P", "Likely pathogenic": "LP", "Uncertain significance": "VUS",
              "Likely benign": "LB", "Benign": "B"}
