import json

from clavagen.cohort import run_cohort
from clavagen.examples import HBB_CDS, write_examples

IU = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M"}


def het(seq, *changes):
    s = list(seq)
    for pos, alt in changes:
        s[pos - 1] = IU["".join(sorted(s[pos - 1] + alt))]
    return "".join(s)


def _run(tmp_path, patients):
    ex = write_examples(str(tmp_path / "in"))
    ann = json.load(open(ex["annotation"]))
    summ = run_cohort(patients, ex["reference_dna"], str(tmp_path / "out"), annotation=ann, online=False)
    folder = tmp_path / "out" / summ["name"]
    return summ, folder


def _vars(folder, pid):
    js = json.load(open(folder / pid / "clinical_summary.json"))
    return {v["hgvs_c"]: v for s in js["samples"] for v in s["variants"]}, js


def test_trio_compound_het_and_de_novo(tmp_path):
    pts = [
        # proband: HbS (c.20A>T) / beta0 (c.118C>T) compound het + de novo c.250G>A
        {"id": "Child", "family": "FAM1", "relationship": "proband", "affected": "affected", "sex": "female",
         "parentage_confirmed": True, "query": het(HBB_CDS, (20, "T"), (118, "T"), (250, "A"))},
        {"id": "Mum", "family": "FAM1", "relationship": "mother", "affected": "unaffected", "sex": "female",
         "query": het(HBB_CDS, (20, "T"))},
        {"id": "Dad", "family": "FAM1", "relationship": "father", "affected": "unaffected", "sex": "male",
         "query": het(HBB_CDS, (118, "T"))},
        {"id": "Sib", "family": "FAM1", "relationship": "sibling", "affected": "unaffected", "sex": "male",
         "query": het(HBB_CDS, (20, "T"))},
    ]
    summ, folder = _run(tmp_path, pts)
    fam = summ["families"][0]
    assert fam["proband"] == "Child" and fam.get("compound_het")
    codes = {(e["code"], e["strength"]) for e in fam["evidence"]}
    assert ("PS2", "Strong") in codes and ("PM3", "Moderate") in codes
    v, js = _vars(folder, "Child")
    assert "PS2" in v["c.250G>A"]["criteria_met"]
    assert v["c.250G>A"]["acmg_classification"] in ("Likely pathogenic", "Pathogenic")   # was VUS alone
    assert "PM3" in v["c.20A>T"]["criteria_met"]
    st = js["samples"][0]["clinical_summary"]["status"]
    assert "Compound heterozygous, in trans (confirmed" in st and "biallelic P/LP" in st
    # parents as a couple: both carriers of recessive HBB P/LP variants -> 25 %
    couple = next(c for c in fam["couples"] if set(c["partners"]) == {"Mum", "Dad"})
    assert couple["risks"][0]["risk"] == "25%"
    # same variant gets the same evidence in the mother
    vm, _ = _vars(folder, "Mum")
    assert "PM3" in vm["c.20A>T"]["criteria_met"]
    assert (folder / "family_FAM1_pedigree.svg").exists()
    html = (folder / "cohort_summary.html").read_text()
    assert "Family and couple analysis" in html and "Segregation table" in html


def test_couple_risk_hbs_hbc(tmp_path):
    pts = [
        {"id": "A", "family": "C1", "relationship": "proband", "sex": "female", "query": het(HBB_CDS, (20, "T"))},
        {"id": "B", "family": "C1", "relationship": "partner", "sex": "male", "query": het(HBB_CDS, (19, "A"))},
    ]
    summ, folder = _run(tmp_path, pts)
    risk = summ["families"][0]["couples"][0]["risks"][0]
    assert risk["risk"] == "25%" and "c.20A>T" in risk["text"] and "c.19G>A" in risk["text"]
    _, js = _vars(folder, "A")
    assert any("Reproductive risk" in r for r in js["samples"][0]["clinical_summary"]["recommendations"])


def test_cis_gives_bp2_and_no_family_without_id(tmp_path):
    pts = [
        {"id": "P", "family": "F", "relationship": "proband", "affected": "affected",
         "query": het(HBB_CDS, (118, "T"), (250, "A"))},
        {"id": "M", "family": "F", "relationship": "mother", "affected": "unaffected",
         "query": het(HBB_CDS, (118, "T"), (250, "A"))},
        {"id": "D", "family": "F", "relationship": "father", "affected": "unaffected", "query": HBB_CDS},
        {"id": "Solo", "query": HBB_CDS},
    ]
    summ, folder = _run(tmp_path, pts)
    assert len(summ["families"]) == 1
    v, _ = _vars(folder, "P")
    assert "BP2" in v["c.250G>A"]["criteria_met"]


def test_de_novo_in_sibling(tmp_path):
    """V11: a variant absent in both parents is reported as de novo in a sibling, not only the proband."""
    pts = [
        {"id": "Pro", "family": "FAM2", "relationship": "proband", "affected": "affected",
         "query": het(HBB_CDS, (20, "T"))},
        {"id": "Mum", "family": "FAM2", "relationship": "mother", "affected": "unaffected",
         "query": het(HBB_CDS, (20, "T"))},
        {"id": "Dad", "family": "FAM2", "relationship": "father", "affected": "unaffected", "query": HBB_CDS},
        {"id": "SibA", "family": "FAM2", "relationship": "sibling", "affected": "affected",
         "query": het(HBB_CDS, (250, "A"))},
        {"id": "SibU", "family": "FAM2", "relationship": "sibling", "affected": "unaffected",
         "query": het(HBB_CDS, (300, "A"))},
    ]
    summ, folder = _run(tmp_path, pts)
    found = " ".join(summ["families"][0]["findings"])
    assert "de novo in sibling SibA" in found and "de novo in sibling SibU" in found
    v, _ = _vars(folder, "SibA")
    assert "PM6" in v["c.250G>A"]["criteria_met"]
    vu, _ = _vars(folder, "SibU")
    assert not any(c in vu[k]["criteria_met"] for k in vu for c in ("PM6", "PS2"))
