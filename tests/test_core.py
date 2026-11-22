import os

from clavagen.acmg import Criterion, combine
from clavagen.align import align
from clavagen.caller import call_variants
from clavagen.consequence import TranscriptModel, annotate_dna
from clavagen.examples import HBB_CDS
from clavagen.genetic_code import BLOSUM62, grantham
from clavagen.pipeline import Pipeline
from clavagen.examples import write_examples


def _call(q):
    a = align("ref", HBB_CDS, "q", q)
    vs, _ = call_variants(a, HBB_CDS)
    tm = TranscriptModel(HBB_CDS, exon_boundaries_c=[92, 315])
    return [annotate_dna(v, tm) for v in vs], vs


def test_blosum_symmetric():
    assert all(BLOSUM62[(a, b)] == BLOSUM62[(b, a)] for a, b in BLOSUM62)


def test_grantham_values():
    assert grantham("S", "R") in (109, 110) and grantham("E", "V") in (121, 122) and grantham("L", "L") == 0


def test_sickle_het():
    q = HBB_CDS[:19] + "W" + HBB_CDS[20:]
    ann, raw = _call(q)
    assert raw[0].zygosity == "heterozygous"
    assert ann[0]["hgvs_c"] == "c.20A>T" and ann[0]["hgvs_p"] == "p.(Glu7Val)"


def test_nonsense_and_nmd():
    ann, _ = _call(HBB_CDS[:117] + "T" + HBB_CDS[118:])
    assert ann[0]["consequence"] == "stop_gained" and ann[0]["hgvs_p"] == "p.(Gln40Ter)"


def test_frameshift_3prime_normalisation():
    ann, _ = _call(HBB_CDS[:123] + HBB_CDS[127:])
    assert ann[0]["hgvs_c"] == "c.126_129del"
    assert ann[0]["hgvs_p"] == "p.(Phe42LeufsTer19)"


def test_inframe_dup():
    q = HBB_CDS[:54] + "AAG" + HBB_CDS[54:]  # duplicate codon 18 (Lys)
    ann, _ = _call(q)
    assert ann[0]["consequence"] == "inframe_insertion" and "dup" in ann[0]["hgvs_p"]


def test_start_lost_and_synonymous():
    ann, _ = _call("ACG" + HBB_CDS[3:])
    assert ann[0]["consequence"] == "start_lost"
    ann, _ = _call(HBB_CDS[:8] + "C" + HBB_CDS[9:])
    assert ann[0]["consequence"] == "synonymous_variant"


def test_combining_rules():
    c = [Criterion("PVS1", "met", "VeryStrong"), Criterion("PM2", "met", "Supporting")]
    assert combine(c)["classification"] == "Likely pathogenic"
    c = [Criterion("BA1", "met", "StandAlone")]
    assert combine(c)["classification"] == "Benign"
    c = [Criterion("PS3", "met", "Strong"), Criterion("BS1", "met", "Strong")]
    assert combine(c)["classification"] == "Uncertain significance"


def test_pipeline_demo(tmp_path):
    p = write_examples(str(tmp_path))
    res = Pipeline([p["dna_queries"], p["protein_query"]], [p["reference_dna"]], annotation=p["annotation"]).run()
    assert len(res["samples"]) == 5
    cls = {v["hgvs_c"]: v["classification"] for s in res["samples"] for v in s["variants"]}
    assert cls["c.118C>T"] in ("Pathogenic", "Likely pathogenic")
    assert cls["c.9T>C"] in ("Likely benign", "Benign")
