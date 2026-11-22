"""Family and couple analysis.

Patients that share a *family ID* are analysed together. Each member has a relationship
to the proband (proband, mother, father, sibling, child, partner, other) and an affected
status. From the members' genotypes the module derives:

* de novo status         -> PS2 (parentage confirmed) / PM6 (assumed)
* phase in the proband   -> PM3 (in trans with a P/LP variant, recessive genes),
                            BP2 (in cis with a P/LP variant)
* segregation            -> PP1 (affected relatives carry it; simplified counting),
                            BS4 (an affected relative does not carry it; dominant genes)
* couple risk            -> reproductive risk for proband+partner and mother+father

Family evidence is applied to the variant in every member who carries it, so the
same variant gets the same classification across the family. Reviewer evidence
always takes precedence.
"""

from __future__ import annotations

from html import escape

from .annotation import normalize_hgvs

ROLES = ["proband", "mother", "father", "sibling", "child", "partner", "other"]
PLP = ("Pathogenic", "Likely pathogenic")
CLASS_COLORS = {"Pathogenic": "#b2182b", "Likely pathogenic": "#ef6548", "Uncertain significance": "#e3a008",
                "Likely benign": "#6baed6", "Benign": "#2166ac"}


def vkey(v):
    """Identity of a variant across family members."""
    name = v.get("hgvs_c") or v.get("hgvs_p") or v.get("hgvs_ref")
    return f"{v.get('ref_id')}|{v.get('mol_type')}|{normalize_hgvs(name)}"


def _inh(ctx_inh):
    s = (ctx_inh or "").upper()
    if "AR" in s and "AD" not in s:
        return "AR"
    if "AD" in s and "AR" not in s:
        return "AD"
    if s.startswith("X") or "XL" in s:
        return "XL"
    if "AD" in s and "AR" in s:
        return "AD/AR"
    return ""


def group_families(patients):
    fams = {}
    for p in patients:
        fid = (p.get("family") or "").strip()
        if fid:
            fams.setdefault(fid, []).append(p)
    return {k: v for k, v in fams.items() if len(v) >= 2}


def _member_variants(res):
    """{vkey: (variant, sample)} plus coverage list [(ref_id, mode, start1, end1)]."""
    vs, cov = {}, []
    for s in res["samples"]:
        a = s["alignment"]
        cov.append((s["reference"].id, s["mode"], a.ref_start + 1, a.ref_end))
        for v in s["variants"]:
            vs[vkey(v)] = (v, s)
    return vs, cov


def _covered(cov, v):
    for ref_id, mode, a, b in cov:
        if ref_id == v.get("ref_id") and mode == v.get("mol_type") and a <= v["ref_pos"] <= b:
            return True
    return False


def _genotype(member_vs, member_cov, key, v):
    if key in member_vs:
        z = member_vs[key][0]["zygosity"]
        return "het" if z == "heterozygous" else ("hom" if z.startswith("homozygous") else "carrier")
    return "absent" if _covered(member_cov, v) else "n/c"


def analyse_family(fid, members, results):
    """members: patient dicts; results: {pid: pipeline result}. Returns analysis dict."""
    ok = [m for m in members if m["id"] in results]
    by_role = {}
    for m in ok:
        by_role.setdefault(m.get("relationship") or "other", []).append(m)
    proband = (by_role.get("proband") or [m for m in ok if m.get("affected") == "affected"] or ok or [None])[0]
    out = {"family": fid, "members": [], "proband": proband["id"] if proband else None, "variants": [],
           "findings": [], "couples": [], "evidence": {}, "warnings": []}
    if not proband:
        out["warnings"].append("No analysable members.")
        return out
    if not by_role.get("proband"):
        out["warnings"].append(f"No member marked as proband - {proband['id']} used as index case.")
    mv = {m["id"]: _member_variants(results[m["id"]]) for m in ok}
    ctx_by_ref = {}
    for m in ok:
        for s in results[m["id"]]["samples"]:
            ctx_by_ref[s["reference"].id] = s["context"]
    for m in ok:
        out["members"].append({"id": m["id"], "relationship": m.get("relationship") or "other",
                               "affected": m.get("affected") or "unknown", "sex": m.get("sex", ""),
                               "parentage_confirmed": bool(m.get("parentage_confirmed"))})
    # union of variants in the family
    allv = {}
    for pid, (vs, _) in mv.items():
        for k, (v, s) in vs.items():
            allv.setdefault(k, (v, s))
    for k, (v, s) in sorted(allv.items(), key=lambda kv: (kv[1][0].get("gene", ""), kv[1][0]["ref_pos"])):
        row = {"key": k, "gene": v.get("gene") or "", "hgvs_c": v.get("hgvs_c") or "", "hgvs_p": v.get("hgvs_p") or "",
               "classification": v["classification"], "genotypes": {}}
        for m in ok:
            vs, cov = mv[m["id"]]
            row["genotypes"][m["id"]] = _genotype(vs, cov, k, v)
        out["variants"].append(row)

    ev = {}  # vkey -> list of (code, strength, rationale)

    def add(key, code, strength, why):
        lst = ev.setdefault(key, [])
        if not any(c == code for c, _, _ in lst):
            lst.append((code, strength, why))

    pvs, _ = mv[proband["id"]]
    mothers, fathers = by_role.get("mother", []), by_role.get("father", [])
    mom, dad = (mothers[0] if mothers else None), (fathers[0] if fathers else None)
    relatives = [m for m in ok if m["id"] != proband["id"]]
    prob_aff = proband.get("affected") != "unaffected"

    for k, (v, s) in pvs.items():
        name = " ".join(x for x in (v.get("hgvs_c"), v.get("hgvs_p")) if x)
        inh = _inh(s["context"].get("inheritance"))
        gm = _genotype(*mv[mom["id"]], k, v) if mom else None
        gf = _genotype(*mv[dad["id"]], k, v) if dad else None
        # ---- de novo
        if mom and dad:
            if gm == "absent" and gf == "absent":
                if prob_aff:
                    if proband.get("parentage_confirmed"):
                        add(k, "PS2", "Strong", f"De novo in {proband['id']}: absent in mother ({mom['id']}) and "
                                                f"father ({dad['id']}); maternity and paternity confirmed.")
                    else:
                        add(k, "PM6", "Moderate", f"Assumed de novo in {proband['id']}: absent in both parents "
                                                  "(parentage not confirmed).")
                out["findings"].append(f"{name} ({v.get('gene')}): de novo in {proband['id']} (absent in both parents).")
            elif gm in ("het", "hom", "carrier") or gf in ("het", "hom", "carrier"):
                src = " and ".join(x for x, g in (("mother", gm), ("father", gf)) if g in ("het", "hom", "carrier"))
                out["findings"].append(f"{name} ({v.get('gene')}): inherited from {src}.")
            elif "n/c" in (gm, gf):
                out["findings"].append(f"{name}: parental origin not determinable (region not covered in a parent).")
        elif mom or dad:
            par, g = (mom, gm) if mom else (dad, gf)
            out["findings"].append(f"{name} ({v.get('gene')}): {'present' if g in ('het', 'hom', 'carrier') else g} "
                                   f"in {par['relationship']} {par['id']} (other parent not tested).")
        # ---- AD: inherited from unaffected parent
        if inh == "AD" and v["classification"] in PLP:
            for par, g in ((mom, gm), (dad, gf)):
                if par and g in ("het", "hom", "carrier") and par.get("affected") == "unaffected":
                    out["findings"].append(f"{name}: P/LP variant in a dominant gene inherited from unaffected "
                                           f"{par['relationship']} {par['id']} - consider reduced penetrance / variable "
                                           "expressivity (BS2 only if the disorder is fully penetrant at that age).")
        # ---- segregation (relatives only)
        aff_car, aff_non, unaff_car = [], [], []
        for m in relatives:
            if m.get("relationship") == "partner":
                continue
            g = _genotype(*mv[m["id"]], k, v)
            if m.get("affected") == "affected":
                (aff_car if g in ("het", "hom", "carrier") else aff_non if g == "absent" else []).append(m["id"])
            elif m.get("affected") == "unaffected" and g in ("het", "hom", "carrier"):
                unaff_car.append(m["id"])
        if prob_aff and inh in ("AD", "XL", "AD/AR") and aff_non:
            add(k, "BS4", "Strong", f"Lack of segregation: affected relative(s) {', '.join(aff_non)} do not carry the "
                                    "variant (consider phenocopy / second cause before relying on BS4).")
            out["findings"].append(f"{name}: does NOT segregate - absent in affected {', '.join(aff_non)} (BS4).")
        elif prob_aff and len(aff_car) >= 2 and not aff_non:
            st = "Strong" if len(aff_car) >= 6 else "Moderate" if len(aff_car) >= 4 else "Supporting"
            add(k, "PP1", st, f"Co-segregates with disease: present in {len(aff_car)} affected relative(s) "
                              f"({', '.join(aff_car)}). Strength from a simplified count - verify informative meioses.")
            out["findings"].append(f"{name}: co-segregates in {len(aff_car)} affected relatives (PP1_{st}).")
        elif aff_car:
            out["findings"].append(f"{name}: also present in affected {', '.join(aff_car)}.")
        if unaff_car and inh == "AD" and v["classification"] in PLP:
            out["findings"].append(f"{name}: present in unaffected relative(s) {', '.join(unaff_car)}.")

    # ---- de novo variants in siblings (same parents as the proband)
    if mom and dad:
        for sib in (m for m in ok if m.get("relationship") == "sibling"):
            svs, _ = mv[sib["id"]]
            for k, (v, s) in svs.items():
                if k in pvs:
                    continue                      # already assessed in the proband
                gm = _genotype(*mv[mom["id"]], k, v)
                gf = _genotype(*mv[dad["id"]], k, v)
                if gm != "absent" or gf != "absent":
                    continue
                name = " ".join(x for x in (v.get("hgvs_c"), v.get("hgvs_p")) if x)
                if sib.get("affected") == "affected":
                    if sib.get("parentage_confirmed"):
                        add(k, "PS2", "Strong", f"De novo in {sib['id']}: absent in mother ({mom['id']}) and father "
                                                f"({dad['id']}); maternity and paternity confirmed.")
                    else:
                        add(k, "PM6", "Moderate", f"Assumed de novo in {sib['id']}: absent in both parents "
                                                  "(parentage not confirmed).")
                    tail = ""
                else:
                    tail = " - sibling is unaffected, so PS2/PM6 are not applied (they need an affected carrier)"
                out["findings"].append(f"{name} ({v.get('gene')}): de novo in sibling {sib['id']} "
                                       f"(absent in both parents){tail}.")

    # ---- phase in the proband (pairs in the same gene)
    pv = [(k, v) for k, (v, s) in pvs.items() if v.get("gene")]
    for i in range(len(pv)):
        for j in range(i + 1, len(pv)):
            (ka, a), (kb, b) = pv[i], pv[j]
            if a["gene"] != b["gene"] or a["zygosity"] != "heterozygous" or b["zygosity"] != "heterozygous":
                continue
            ga = (_genotype(*mv[mom["id"]], ka, a) if mom else None, _genotype(*mv[dad["id"]], ka, a) if dad else None)
            gb = (_genotype(*mv[mom["id"]], kb, b) if mom else None, _genotype(*mv[dad["id"]], kb, b) if dad else None)
            car = lambda g: g in ("het", "hom", "carrier")  # noqa: E731
            phase, how = None, ""
            if mom and dad:
                if (car(ga[0]) and not car(ga[1]) and car(gb[1]) and not car(gb[0])) or \
                        (car(ga[1]) and not car(ga[0]) and car(gb[0]) and not car(gb[1])):
                    phase, how = "trans", "confirmed by testing both parents"
                elif (car(ga[0]) and car(gb[0]) and not car(ga[1]) and not car(gb[1])) or \
                        (car(ga[1]) and car(gb[1]) and not car(ga[0]) and not car(gb[0])):
                    phase, how = "cis", "both inherited from the same parent"
            elif mom or dad:
                idx = 0 if mom else 1
                if car(ga[idx]) != car(gb[idx]) and "n/c" not in (ga[idx], gb[idx]):
                    phase, how = "trans", "inferred from one parent only"
                elif car(ga[idx]) and car(gb[idx]):
                    phase, how = "cis", "both present in the tested parent"
            if not phase:
                continue
            na = a.get("hgvs_c") or a.get("hgvs_p")
            nb = b.get("hgvs_c") or b.get("hgvs_p")
            out["findings"].append(f"{a['gene']}: {na} and {nb} are in {phase} in {proband['id']} ({how}).")
            out.setdefault("phases", []).append({"gene": a["gene"], "a": ka, "b": kb, "phase": phase, "how": how,
                                                 "member": proband["id"]})
            inh = _inh(ctx_by_ref.get(a.get("ref_id"), {}).get("inheritance"))
            for (kx, x), (ky, y) in (((ka, a), (kb, b)), ((kb, b), (ka, a))):
                if y["classification"] not in PLP:
                    continue
                ny = y.get("hgvs_c") or y.get("hgvs_p")
                if phase == "trans" and inh in ("AR", "AD/AR"):
                    st = "Moderate" if how.startswith("confirmed") else "Supporting"
                    add(kx, "PM3", st, f"In trans with {y['classification'].lower()} variant {ny} in "
                                       f"{proband['id']} ({how}).")
                elif phase == "cis":
                    add(kx, "BP2", "Supporting", f"In cis with {y['classification'].lower()} variant {ny} ({how}).")
            if phase == "trans" and a["classification"] in PLP and b["classification"] in PLP and inh in ("AR", "AD/AR"):
                out["findings"].append(f"Compound heterozygous P/LP variants in {a['gene']} in {proband['id']} "
                                       "(phase established) - consistent with the recessive condition.")
                out["compound_het"] = True

    # ---- couples
    couples = []
    for p in by_role.get("partner", []):
        couples.append((proband, p))
    if mom and dad:
        couples.append((mom, dad))
    for x, y in couples:
        out["couples"].append(couple_risk(x, y, results[x["id"]], results[y["id"]]))
    out["evidence"] = ev
    return out


def _plp_by_gene(res):
    g = {}
    for s in res["samples"]:
        for v in s["variants"]:
            if v["classification"] in PLP and v.get("gene"):
                g.setdefault(v["gene"], []).append((v, _inh(s["context"].get("inheritance")),
                                                    s["context"].get("diseases") or []))
    return g


def couple_risk(a, b, ra, rb):
    """Reproductive risk for a couple from their P/LP variants."""
    ga, gb = _plp_by_gene(ra), _plp_by_gene(rb)
    risks = []
    for gene in sorted(set(ga) | set(gb)):
        va, vb = ga.get(gene, []), gb.get(gene, [])
        inh = (va or vb)[0][1]
        dis = ", ".join(d["name"] for d in (va or vb)[0][2][:2]) or f"{gene}-related disorder"
        name = lambda lst: " / ".join((v.get("hgvs_c") or v.get("hgvs_p")) for v, _, _ in lst)  # noqa: E731
        hom = lambda lst: any(v["zygosity"].startswith("homozygous") for v, _, _ in lst)  # noqa: E731
        if inh in ("AR", "AD/AR") and va and vb:
            pa, pb = (1.0 if hom(va) else 0.5), (1.0 if hom(vb) else 0.5)
            risk = pa * pb
            risks.append({"gene": gene, "level": "high", "risk": f"{risk * 100:.0f}%",
                          "text": f"Both partners carry P/LP variants in {gene} ({a['id']}: {name(va)}; "
                                  f"{b['id']}: {name(vb)}). Each pregnancy has a {risk * 100:.0f}% risk of an affected "
                                  f"child with {dis} (biallelic: {name(va)} + {name(vb)}). Offer genetic counselling "
                                  "and reproductive options (prenatal / preimplantation testing)."})
        elif inh in ("AR", "AD/AR") and (va or vb):
            who, lst = (a, va) if va else (b, vb)
            other = b if va else a
            risks.append({"gene": gene, "level": "low", "risk": "low",
                          "text": f"Only {who['id']} carries a P/LP variant in recessive {gene} ({name(lst)}). Risk of an "
                                  f"affected child is low if {other['id']}'s analysed region is complete; children have "
                                  "a 50% chance of being carriers. Note: a partner's negative result covers only the "
                                  "sequenced region."})
        elif inh == "AD":
            for who, lst in ((a, va), (b, vb)):
                if lst:
                    risks.append({"gene": gene, "level": "high", "risk": "50%",
                                  "text": f"{who['id']} carries a P/LP variant in dominant {gene} ({name(lst)}): each child "
                                          f"has a 50% chance of inheriting it ({dis}; penetrance may be incomplete)."})
        elif inh == "XL":
            for who, lst in ((a, va), (b, vb)):
                if not lst:
                    continue
                if (who.get("sex") or "").lower() == "female":
                    risks.append({"gene": gene, "level": "high", "risk": "50% of sons",
                                  "text": f"{who['id']} (female) carries an X-linked P/LP variant in {gene} ({name(lst)}): "
                                          "each son has a 50% risk of being affected and each daughter a 50% chance of "
                                          "being a carrier."})
                else:
                    risks.append({"gene": gene, "level": "high", "risk": "all daughters carriers",
                                  "text": f"{who['id']} carries an X-linked P/LP variant in {gene} ({name(lst)}): all "
                                          "daughters will be carriers; sons are not affected through him."})
    if not risks:
        risks.append({"gene": "", "level": "none", "risk": "-",
                      "text": "No P/LP variants detected in either partner in the analysed regions."})
    return {"partners": [a["id"], b["id"]], "label": f"{a['id']} + {b['id']}", "risks": risks}


# ------------------------------------------------------------------ pedigree drawing
def pedigree_svg(fa, width=760):
    """Simple 3-generation pedigree with genotypes for the family's variants of interest."""
    mem = {m["id"]: m for m in fa["members"]}
    proband = fa["proband"]
    gen1 = [m for m in fa["members"] if m["relationship"] in ("mother", "father")]
    gen1.sort(key=lambda m: 0 if m["relationship"] == "father" else 1)
    gen2 = [mem[proband]] + [m for m in fa["members"] if m["relationship"] in ("sibling", "partner")]
    gen3 = [m for m in fa["members"] if m["relationship"] == "child"]
    others = [m for m in fa["members"] if m["relationship"] == "other"]
    show = [v for v in fa["variants"] if v["classification"] in PLP or v["classification"] == "Uncertain significance"][:3]
    if not show:
        show = fa["variants"][:3]
    rows = [r for r in (gen1, gen2, gen3, others) if r]
    gh = 150 + 14 * len(show)
    H = 40 + gh * len(rows) + 16
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{H}" viewBox="0 0 {width} {H}">'
             '<rect width="100%" height="100%" fill="#fff"/>',
             f'<text x="12" y="22" font-size="14" font-weight="bold" font-family="Helvetica,Arial">Family {escape(fa["family"])}'
             '</text>']
    pos = {}
    S = 34

    def sym(m, x, y):
        aff = m["affected"] == "affected"
        fill = "#13294b" if aff else ("#fff" if m["affected"] == "unaffected" else "url(#unk)")
        sex = (m.get("sex") or "").lower()
        if sex == "male":
            shape = f'<rect x="{x - S / 2}" y="{y - S / 2}" width="{S}" height="{S}" fill="{fill}" stroke="#13294b" stroke-width="2"/>'
        elif sex == "female":
            shape = f'<circle cx="{x}" cy="{y}" r="{S / 2}" fill="{fill}" stroke="#13294b" stroke-width="2"/>'
        else:
            h = S / 2 + 3
            shape = (f'<polygon points="{x},{y - h} {x + h},{y} {x},{y + h} {x - h},{y}" fill="{fill}" '
                     'stroke="#13294b" stroke-width="2"/>')
        s = shape
        if m["id"] == proband:
            # proband arrow, lifted so it points at the symbol and clears the name label below it
            s += (f'<g transform="translate(-4,{-S / 2})">'
                  f'<line x1="{x - S}" y1="{y + S}" x2="{x - S / 2 - 3}" y2="{y + S / 2 + 3}" stroke="#b2182b" stroke-width="2.5"/>'
                  f'<polygon points="{x - S / 2 - 3},{y + S / 2 + 3} {x - S / 2 - 12},{y + S / 2 + 6} {x - S / 2 - 6},{y + S / 2 + 12}" fill="#b2182b"/></g>')
        s += (f'<text x="{x}" y="{y + S / 2 + 16}" font-size="11.5" font-weight="bold" text-anchor="middle" '
              f'font-family="Helvetica,Arial">{escape(m["id"])}</text>'
              f'<text x="{x}" y="{y + S / 2 + 29}" font-size="10" text-anchor="middle" fill="#666" '
              f'font-family="Helvetica,Arial">{escape(m["relationship"])}</text>')
        for i, v in enumerate(show):
            g = next((r["genotypes"].get(m["id"]) for r in fa["variants"] if r["key"] == v["key"]), "n/c")
            lab = {"het": "+/-", "hom": "+/+", "absent": "-/-", "carrier": "+", "n/c": "n.t."}.get(g, g)
            col = CLASS_COLORS.get(v["classification"], "#444") if g in ("het", "hom", "carrier") else "#999"
            s += (f'<text x="{x}" y="{y + S / 2 + 44 + i * 14}" font-size="10.5" text-anchor="middle" fill="{col}" '
                  f'font-family="Menlo,Consolas,monospace">{escape(lab)}</text>')
        return s

    for gi, row in enumerate(rows):
        y = 40 + gi * gh + 40
        n = len(row)
        for i, m in enumerate(row):
            x = width * (i + 1) / (n + 1)
            pos[m["id"]] = (x, y)
    # relationship lines
    lines = []
    fa_m = [m for m in gen1]
    if len(fa_m) == 2:
        (x1, y1), (x2, y2) = pos[fa_m[0]["id"]], pos[fa_m[1]["id"]]
        lines.append(f'<line x1="{x1 + S / 2}" y1="{y1}" x2="{x2 - S / 2}" y2="{y2}" stroke="#13294b" stroke-width="2"/>')
        mx = (x1 + x2) / 2
        kids = [m for m in gen2 if m["relationship"] in ("sibling",) or m["id"] == proband]
        if kids:
            ky = pos[kids[0]["id"]][1] - S / 2 - 20
            lines.append(f'<line x1="{mx}" y1="{y1}" x2="{mx}" y2="{ky}" stroke="#13294b" stroke-width="2"/>')
            xs = [pos[k["id"]][0] for k in kids]
            lines.append(f'<line x1="{min(xs + [mx])}" y1="{ky}" x2="{max(xs + [mx])}" y2="{ky}" stroke="#13294b" stroke-width="2"/>')
            for k in kids:
                kx, kyy = pos[k["id"]]
                lines.append(f'<line x1="{kx}" y1="{ky}" x2="{kx}" y2="{kyy - S / 2 - 2}" stroke="#13294b" stroke-width="2"/>')
    elif len(fa_m) == 1:
        (x1, y1) = pos[fa_m[0]["id"]]
        px, py = pos[proband]
        lines.append(f'<line x1="{x1}" y1="{y1 + S / 2 + 48 + 14 * len(show)}" x2="{px}" y2="{py - S / 2 - 2}" '
                     'stroke="#13294b" stroke-width="2" stroke-dasharray="4,3"/>')
    partners = [m for m in gen2 if m["relationship"] == "partner"]
    if partners:
        px, py = pos[proband]
        qx, qy = pos[partners[0]["id"]]
        lines.append(f'<line x1="{min(px, qx) + S / 2}" y1="{py}" x2="{max(px, qx) - S / 2}" y2="{qy}" stroke="#13294b" stroke-width="2"/>')
        if gen3:
            mx = (px + qx) / 2
            ky = pos[gen3[0]["id"]][1] - S / 2 - 20
            lines.append(f'<line x1="{mx}" y1="{py}" x2="{mx}" y2="{ky}" stroke="#13294b" stroke-width="2"/>')
            xs = [pos[c["id"]][0] for c in gen3]
            lines.append(f'<line x1="{min(xs + [mx])}" y1="{ky}" x2="{max(xs + [mx])}" y2="{ky}" stroke="#13294b" stroke-width="2"/>')
            for c in gen3:
                cx, cy = pos[c["id"]]
                lines.append(f'<line x1="{cx}" y1="{ky}" x2="{cx}" y2="{cy - S / 2 - 2}" stroke="#13294b" stroke-width="2"/>')
    elif gen3:
        px, py = pos[proband]
        for c in gen3:
            cx, cy = pos[c["id"]]
            lines.append(f'<line x1="{px}" y1="{py + S / 2 + 48 + 14 * len(show)}" x2="{cx}" y2="{cy - S / 2 - 2}" '
                         'stroke="#13294b" stroke-width="2"/>')
    parts.append('<defs><pattern id="unk" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
                 '<line x1="0" y1="0" x2="0" y2="6" stroke="#9aa5b8" stroke-width="2"/></pattern></defs>')
    parts += lines
    for m in fa["members"]:
        if m["id"] in pos:
            parts.append(sym(m, *pos[m["id"]]))
    leg_y = H - 30
    legend = "Genotype rows: " + "; ".join(f"{i + 1}) {v['gene']} {v['hgvs_c'] or v['hgvs_p']} [{v['classification']}]"
                                           for i, v in enumerate(show))
    parts.append(f'<text x="12" y="{leg_y - 14}" font-size="10" fill="#444" font-family="Helvetica,Arial">'
                 f'{escape(legend[:150])}</text>')
    parts.append(f'<text x="12" y="{leg_y - 14 + 14}" font-size="10" fill="#666" font-family="Helvetica,Arial">'
                 'Filled = affected, open = unaffected, hatched = unknown; square = male, circle = female, diamond = '
                 'unknown sex; arrow = proband.</text>')
    parts.append(f'<text x="12" y="{leg_y + 14}" font-size="10" fill="#666" font-family="Helvetica,Arial">'
                 'Genotypes: +/- heterozygous, +/+ homozygous, -/- absent, n.t. = not tested / region not covered.</text>')
    parts.append("</svg>")
    return "".join(parts)


def family_html(fa):
    """HTML block: findings, couple risks, segregation table and pedigree."""
    mem_ids = [m["id"] for m in fa["members"]]
    head = "".join(f"<th>{escape(i)}<br><span class='muted'>{escape(next(m['relationship'] for m in fa['members'] if m['id'] == i))}"
                   f", {escape(next(m['affected'] for m in fa['members'] if m['id'] == i))}</span></th>" for i in mem_ids)
    sym = {"het": "&#9680; het", "hom": "&#9679; hom", "absent": "&#9675; absent", "carrier": "&#9679; carrier", "n/c": "n.t."}
    rows = "".join(
        f"<tr><td><b>{escape(v['gene'])}</b> <span class='mono'>{escape(v['hgvs_c'])} {escape(v['hgvs_p'])}</span><br>"
        f"<span class='badge' style='background:{CLASS_COLORS.get(v['classification'], '#999')}'>{escape(v['classification'])}</span></td>"
        + "".join(f"<td>{sym.get(v['genotypes'].get(i, 'n/c'), escape(v['genotypes'].get(i, '')))}</td>" for i in mem_ids)
        + "</tr>" for v in fa["variants"])
    couples = ""
    for c in fa["couples"]:
        items = "".join(f"<li class='{ 'risk' if r['level'] == 'high' else ''}'><b>{escape(r['risk'])}</b> - "
                        f"{escape(r['text'])}</li>" for r in c["risks"])
        couples += f"<h4>Couple: {escape(c['label'])}</h4><ul>{items}</ul>"
    ev = ""
    for k, lst in fa["evidence"].items():
        name = next((v["hgvs_c"] or v["hgvs_p"] for v in fa["variants"] if v["key"] == k), k)
        ev += "".join(f"<li><b>{escape(code)}_{escape(st)}</b> for {escape(name)}: {escape(why)}</li>" for code, st, why in lst)
    return (f"<div class='family'><h3>Family {escape(fa['family'])} (proband: {escape(fa['proband'] or '-')})</h3>"
            + "".join(f"<p class='warn'>{escape(w)}</p>" for w in fa["warnings"])
            + f"<div class='ped'>{pedigree_svg(fa)}</div>"
            + "<h4>Findings</h4><ul>" + ("".join(f"<li>{escape(x)}</li>" for x in fa["findings"]) or "<li>None.</li>") + "</ul>"
            + ("<h4>ACMG evidence added from family data</h4><ul>" + ev + "</ul>" if ev else "")
            + couples
            + "<h4>Segregation table</h4><table><tr><th>Variant</th>" + head + "</tr>" + rows + "</table></div>")
