"""V11: frameshifts, premature stops and extensions in protein-level queries."""
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from clavagen.cohort import run_cohort  # noqa: E402
from clavagen.examples import HBB_CDS  # noqa: E402
from clavagen.genetic_code import translate  # noqa: E402

PROT = translate(HBB_CDS).rstrip("*")
FS = translate(HBB_CDS[:125] + HBB_CDS[129:], to_stop=True)          # c.126_129del
INS = translate(HBB_CDS[:50] + "A" + HBB_CDS[50:], to_stop=True)       # c.50_51insA


def _run(tmp_path, ref_seq, ref_name, cases):
    ref = tmp_path / "ref.fa"
    ref.write_text(f">{ref_name}\n{ref_seq}\n")
    pts = [{"id": k, "query": f">{k}\n{v}"} for k, v in cases.items()]
    run_cohort(pts, str(ref), str(tmp_path / "out"), gene="HBB", online=False, log=lambda *a: None)
    res = {}
    for k in cases:
        f = glob.glob(str(tmp_path / "out" / "*" / k / "variants.tsv"))[0]
        res[k] = list(csv.DictReader(open(f), delimiter="\t"))
    return res


CASES = {"fs_stop": FS, "fs_open": FS.rstrip("*"), "ins": INS, "nonsense": PROT[:39] + "*",
         "partial": PROT[:80], "ext": PROT + "GGKL", "missense": PROT[:6] + "V" + PROT[7:]}


def _check(res, dna_ref):
    fs = res["fs_stop"]
    assert len(fs) == 1 and fs[0]["consequence"] == "frameshift_variant"
    assert fs[0]["hgvs_p"] == "p.(Phe42LeufsTer19)"
    assert "PVS1" in fs[0]["criteria_met"]
    assert res["fs_open"][0]["hgvs_p"] == "p.(Phe42LeufsTer?)"
    ins = res["ins"]
    assert len(ins) == 1, "chance-aligned residue must be absorbed into the frameshift"
    assert ins[0]["hgvs_p"] == "p.(Lys18GlnfsTer6)"
    assert res["nonsense"][0]["hgvs_p"] == "p.(Gln40Ter)"
    assert res["nonsense"][0]["consequence"] == "stop_gained"
    assert res["partial"] == []
    assert res["ext"][0]["consequence"] == "stop_lost"
    assert res["missense"][0]["hgvs_p"] == "p.(Glu7Val)" and len(res["missense"]) == 1
    if dna_ref:
        assert fs[0]["hgvs_c_full"].endswith("c.126_129del")
        assert res["nonsense"][0]["hgvs_c_full"].endswith("c.118C>T")
        assert ins[0]["hgvs_c_full"] == ""          # ambiguous: several insertions give this protein


def test_protein_query_dna_reference(tmp_path):
    _check(_run(tmp_path, HBB_CDS, "HBB_NM_000518.5", CASES), dna_ref=True)


def test_protein_query_protein_reference(tmp_path):
    _check(_run(tmp_path, PROT, "HBB_NP_000509.1", CASES), dna_ref=False)


def test_partial_query_warns(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(f">HBB_NP_000509.1\n{PROT}\n")
    run_cohort([{"id": "p", "query": ">p\n" + PROT[:80]}], str(ref), str(tmp_path / "o"), gene="HBB",
               online=False, log=lambda *a: None)
    md = open(glob.glob(str(tmp_path / "o" / "*" / "p" / "clinical_summary.md"))[0]).read()
    assert "ends at residue 80 of 147" in md


def test_frameshift_dna_candidates():
    from clavagen.truncation import infer_frameshift_dna
    v = {"aa_pos": 18, "alt_tail": INS[17:]}
    c = infer_frameshift_dna(HBB_CDS, "", 0, v)
    assert "c.50_51insA" in c and "c.50dup" in c
    v = {"aa_pos": 42, "alt_tail": FS[41:]}
    assert infer_frameshift_dna(HBB_CDS, "", 0, v) == ["c.126_129del"]
