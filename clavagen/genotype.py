"""Patient-level genotype per gene (V12).

A patient's variants can come from several files (amplicons, forward/reverse traces, a consensus).
Interpretation per file misses genotypes that need two variants together, e.g. a compound
heterozygote whose two variants sit in different amplicons. This module merges all of a patient's
variants per gene and describes the genotype:

* zygosity across files (a heterozygous call from a trace beats an "apparent homozygous" consensus);
* phase from the family analysis (in trans / in cis), or "phase unknown";
* the variant type of each allele (null, missense, in-frame, splice region ...);
* a genotype statement that does not depend on the variants being P/LP - two VUS in trans in an
  affected child of a recessive family are reported as a compound heterozygous genotype, with the
  classification of each allele stated next to it.
"""
from .consequence import CONSEQUENCE_LABEL
from .family import vkey

PLP = ("Pathogenic", "Likely pathogenic")
BENIGN = ("Benign", "Likely benign")

NULL = {"stop_gained": "nonsense", "frameshift_variant": "frameshift", "splice_donor_variant": "canonical splice",
        "splice_acceptor_variant": "canonical splice", "start_lost": "start-loss", "exon_loss_variant": "exon loss"}


def variant_type(v):
    """(group, label): group is null | missense | inframe | stoploss | splice_region | synonymous | other."""
    c = v.get("consequence") or ""
    if c in NULL:
        return "null", f"null ({NULL[c]})"
    if c == "missense_variant":
        return "missense", "missense"
    if c in ("inframe_insertion", "inframe_deletion", "protein_altering_variant"):
        return "inframe", CONSEQUENCE_LABEL.get(c, c).lower()
    if c == "stop_lost":
        return "stoploss", "stop-loss (extension)"
    if c == "splice_region_variant":
        return "splice_region", "splice region"
    if c in ("synonymous_variant", "stop_retained_variant"):
        return "synonymous", "synonymous"
    return "other", CONSEQUENCE_LABEL.get(c, c or "other").lower()


def _name(v):
    return " ".join(x for x in (v.get("hgvs_c"), v.get("hgvs_p")) if x) or v.get("hgvs_ref") or "?"


def _merge_zyg(a, b):
    if a == b:
        return a
    if "heterozygous" in (a, b):          # trace/IUPAC heterozygous call wins over an apparent-homozygous consensus
        return "heterozygous"
    return a if not a.startswith("not determinable") else b


def gene_key(sample):
    """Group by gene symbol; when no symbol is known, by reference sequence."""
    return sample["context"].get("gene") or f"[{sample['reference'].id}]"


def merged_variants(samples, gene):
    """De-duplicated variants of one gene across all of a patient's samples."""
    out = {}
    for s in samples:
        if gene_key(s) != gene:
            continue
        for v in s["variants"]:
            k = vkey(v)
            trace = bool(getattr(s["query"], "trace", None))
            if k in out:
                out[k]["zygosity"] = _merge_zyg(out[k]["zygosity"], v["zygosity"])
                out[k]["_queries"].add(s["query"].id)
                out[k]["_from_trace"] = out[k]["_from_trace"] or trace
            else:
                d = dict(v)
                d["_queries"] = {s["query"].id}
                d["_from_trace"] = trace
                out[k] = d
    return out


def is_indel(v):
    name = (v.get("hgvs_c") or "") + " " + (v.get("kind") or "")
    return any(t in name for t in ("del", "ins", "dup", "frameshift", "insertion", "deletion", "duplication"))


def _class_short(c):
    return {"Pathogenic": "P", "Likely pathogenic": "LP", "Uncertain significance": "VUS",
            "Likely benign": "LB", "Benign": "B"}.get(c, c)


def _type_comment(types, gene, lof):
    groups = sorted(t[0] for t in types)
    lines = []
    if groups == ["null", "null"]:
        lines.append("Both alleles are null (loss-of-function) variants: no functional protein is expected from "
                     "either allele - typically the more severe end of the spectrum when loss of function is the "
                     "disease mechanism.")
    elif "null" in groups and len(groups) == 2:
        other = [t[1] for t in types if t[0] != "null"][0]
        lines.append(f"One null (loss-of-function) allele and one {other} allele: the phenotype depends on the "
                     "residual function of the non-null allele, which may be hypomorphic - presentation can be "
                     "milder or more variable than with two null alleles.")
    elif groups == ["missense", "missense"]:
        lines.append("Both alleles are missense variants: some residual protein function is possible; severity is "
                     "difficult to predict from the genotype alone.")
    elif len(groups) == 2:
        lines.append(f"Allele types: {types[0][1]} / {types[1][1]}.")
    if "null" in groups and lof is not True:
        lines.append(f"Loss of function is not established as a disease mechanism for {gene} in the annotation "
                     "(lof_mechanism unknown), so PVS1 could not be applied to the null allele. If LOF is an "
                     "accepted mechanism, set it in the annotation file - the null allele may then reach P/LP.")
    return lines


def gene_genotype(samples, gene, ctx, sex="", phases=None, patient_id=""):
    """Genotype description for one gene. Returns None when no relevant (non-benign) variant is present."""
    allv = merged_variants(samples, gene)
    rel = {k: v for k, v in allv.items() if v["classification"] not in BENIGN}
    if not rel:
        return None
    inh = (ctx.get("inheritance") or "").upper()
    ar = "AR" in inh
    lof = ctx.get("lof_mechanism")
    items = []
    for k, v in sorted(rel.items(), key=lambda kv: kv[1].get("ref_pos", 0)):
        grp, lab = variant_type(v)
        zyg = v["zygosity"]
        # A consensus/text sequence can only show ONE allele of an indel: a heterozygous indel looks
        # homozygous. Without a trace (double peaks) the zygosity of such an indel is not known.
        if zyg.startswith("homozygous") and is_indel(v) and not v.get("_from_trace"):
            zyg = "unresolved (indel in a consensus sequence)"
        items.append({"key": k, "name": _name(v), "type": lab, "type_group": grp,
                      "classification": v["classification"], "zygosity": zyg,
                      "files": sorted(v["_queries"])})
    het = [x for x in items if x["zygosity"] == "heterozygous"]
    hom = [x for x in items if x["zygosity"].startswith("homozygous")]
    nd = [x for x in items if x["zygosity"].startswith("not determinable")]
    unres = [x for x in items if x["zygosity"].startswith("unresolved")]
    het_like = het + unres            # an unresolved indel is at least on one allele
    desc = lambda x: f"{x['name']} - {x['type']}, {x['classification']}"  # noqa: E731
    g = {"gene": gene, "inheritance": inh or "not specified", "variants": items, "phase": None,
         "phase_how": "", "genotype": "", "allele_types": "", "status": None, "lines": [], "recs": []}

    # phase from the family analysis (pairs in this gene for this patient)
    ph = None
    if len(het_like) >= 2 and phases:
        keys = {x["key"] for x in het_like}
        for p in phases:
            if p.get("member") in (patient_id, None) and {p["a"], p["b"]} <= keys:
                if ph is None or p["phase"] == "trans":
                    ph = p
    nplp = lambda xs: sum(1 for x in xs if x["classification"] in PLP)  # noqa: E731

    if hom and (ar or not inh):
        x = hom[0]
        male_x = inh.startswith("X") and sex.lower() in ("male", "m")
        g["genotype"] = "hemizygous" if male_x else "homozygous (apparent)"
        g["allele_types"] = f"{x['type']} / {x['type']}"
        g["lines"].append(f"{gene}: {g['genotype']} {desc(x)}.")
        biallelic = [x, x]
    elif len(het_like) >= 2:
        a, b = het_like[0], het_like[1]
        if ph and ph["phase"] == "trans":
            g["genotype"], g["phase"], g["phase_how"] = "compound heterozygous", "trans", ph.get("how", "")
            g["lines"].append(f"{gene}: compound heterozygous - {desc(a)}; {desc(b)} - in trans ({g['phase_how']}).")
            biallelic = [a, b]
        elif ph and ph["phase"] == "cis":
            g["genotype"], g["phase"], g["phase_how"] = "two heterozygous variants in cis", "cis", ph.get("how", "")
            g["lines"].append(f"{gene}: {desc(a)} and {desc(b)} are in cis (same allele; {g['phase_how']}) - "
                              "the other allele carries neither, so this is not a biallelic genotype.")
            biallelic = None
        else:
            g["genotype"], g["phase"] = "two heterozygous variants, phase unknown (possible compound heterozygote)", "unknown"
            g["lines"].append(f"{gene}: two heterozygous variants - {desc(a)}; {desc(b)}. Phase is not known: "
                              "compound heterozygous if they are on different alleles (in trans).")
            g["recs"].append("Establish phase (parental testing or long-read sequencing) to confirm compound "
                             "heterozygosity.")
            biallelic = None
        g["allele_types"] = f"{a['type']} / {b['type']}"
        if unres:
            g["lines"].append("The zygosity of " + ", ".join(x["name"] for x in unres) + " is unresolved (indel in a "
                              "consensus sequence); confirm with the trace.")
        if len(het_like) > 2:
            g["lines"].append(f"{len(het_like)} heterozygous variants in {gene} in total - phase of the others not assessed.")
    elif len(nd) >= 2:
        g["genotype"] = "two variants, zygosity not determinable (protein input)"
        g["allele_types"] = f"{nd[0]['type']} / {nd[1]['type']}"
        sep = len({f for x in nd for f in x["files"]}) > 1
        g["lines"].append(f"{gene}: {desc(nd[0])}; {desc(nd[1])} (protein input)."
                          + (" They come from separate input sequences: if these represent the patient's two "
                             "alleles, the genotype is compound heterozygous - confirm at DNA level with phase."
                             if sep else ""))
        g["status"] = (f"Two variants (protein input){' on separate sequences - possible compound heterozygote' if sep else ''}"
                       f" [{g['allele_types']}]")
        g["recs"].append("Confirm both variants and their phase at DNA level (zygosity cannot be determined from protein).")
        biallelic = None
    elif unres and not het and not hom:
        x = unres[0]
        g["genotype"] = "indel, zygosity unresolved"
        g["allele_types"] = x["type"]
        g["lines"].append(f"{gene}: {desc(x)} seen in a consensus/text sequence, which can show only one allele - "
                          "a heterozygous indel looks homozygous there. Zygosity is not known.")
        g["recs"].append("Determine the zygosity of the indel from the Sanger trace (.ab1: double peaks after the "
                         "indel) or by allele-specific testing.")
        g["status"] = f"{x['type'].capitalize()} {_class_short(x['classification'])}, zygosity unresolved [{x['type']}]"
        return g
    else:
        x = items[0]
        g["genotype"] = ("heterozygous" if x["zygosity"] == "heterozygous" else
                         "zygosity not determinable" if x["zygosity"].startswith("not") else x["zygosity"])
        g["allele_types"] = x["type"]
        g["lines"].append(f"{gene}: {g['genotype']} {desc(x)}.")
        if x["classification"] == "Uncertain significance":
            g["status_if_none"] = f"{g['genotype'].capitalize()} VUS [{x['type']}]"
        if x["type_group"] == "null" and lof is not True:
            g["lines"] += [
                f"Null variant, but loss of function is not established as a disease mechanism for {gene} in the "
                "annotation (lof_mechanism unknown), so PVS1 was not applied."]
        return g                                       # single variant: the classic interpretation applies

    if biallelic:
        n = nplp(biallelic) if biallelic[0] is not biallelic[1] else (2 if biallelic[0]["classification"] in PLP else 0)
        diag = ", ".join(d["name"] for d in ctx.get("diseases") or []) or f"the {gene}-related condition"
        label = g["genotype"].capitalize() + (f", in trans ({g['phase_how']})" if g["phase"] == "trans" else "")
        types = f"[{g['allele_types']}]"
        if ar:
            if n == 2:
                g["status"] = f"{label} - biallelic P/LP {types}"
                g["lines"].append(f"Biallelic pathogenic genotype - consistent with {diag} if the phenotype is "
                                  "concordant.")
            elif n == 1:
                g["status"] = f"{label} - one P/LP + one VUS {types}"
                g["lines"].append("One allele is P/LP and the other is a VUS: a recessive diagnosis is possible but "
                                  "not established until the VUS is resolved (PM3 supports the VUS when the P/LP "
                                  "variant is in trans).")
                g["recs"].append("Resolve the VUS: functional studies, phenotype specificity, segregation in further "
                                 "relatives, ClinVar/literature re-check.")
            else:
                g["status"] = f"{label} - both variants VUS {types}"
                g["lines"].append(f"The genotype fits autosomal recessive inheritance (biallelic, "
                                  f"{'segregates with the family' if g['phase'] == 'trans' else 'apparent homozygous'})"
                                  ", but both variants are of uncertain significance - not diagnostic on its own.")
                g["recs"].append("Gather evidence for both VUS (functional assays, phenotype match, further "
                                 "segregation, reanalysis); re-run with online annotation for ClinVar/gnomAD data.")
        elif "AD" in inh:
            g["status"] = f"{label} - dominant gene {types}"
            g["lines"].append(f"{gene} is annotated as dominant: one P/LP allele is sufficient for disease; "
                              "biallelic variants may be associated with a different or more severe phenotype "
                              "(gene-specific).")
        else:
            g["status"] = f"{label} {types}"
            g["lines"].append("Inheritance is not specified - set it per reference to interpret the genotype.")
        g["lines"] += _type_comment([(x["type_group"], x["type"]) for x in biallelic], gene, lof)
    elif g["phase"] == "unknown" and ar:
        g["status"] = f"Two heterozygous variants, phase unknown ({g['allele_types']})"
    return g


def apply_patient_genotype(samples, patient, phases=None):
    """Re-interpret each sample using all of the patient's variants in that gene (and family phase)."""
    from .pipeline import clinical_interpretation
    sex = (patient or {}).get("sex") or ""
    pid = (patient or {}).get("id") or ""
    for s in samples:
        gene = gene_key(s)
        merged = list(merged_variants(samples, gene).values())
        fam = s.get("clinical", {}).get("family")
        c = clinical_interpretation(s, patient or {}, variants=merged)
        g = gene_genotype(samples, gene, s["context"], sex, phases, pid)
        if g:
            c["genotype"] = {k: v for k, v in g.items() if k not in ("lines", "recs", "status_if_none")}
            if g["status"]:
                c["status"] = g["status"]
            elif g.get("status_if_none") and c["status"] == "No reportable finding":
                c["status"] = g["status_if_none"]
            c["interpretation"] = g["lines"] + [x for x in c["interpretation"]
                                                 if not x.startswith("Two or more heterozygous P/LP")]
            c["recommendations"] = list(dict.fromkeys(g["recs"] + c["recommendations"]))
        if fam:
            c["family"] = fam
        s["clinical"] = c
    return samples
