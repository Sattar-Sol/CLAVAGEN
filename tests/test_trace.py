import json
import os

from clavagen.examples import HBB_CDS, write_examples
from clavagen.genetic_code import reverse_complement
from clavagen.pipeline import Pipeline
from clavagen.trace import load_trace, read_abif, write_synthetic_ab1

D = os.path.join(os.path.dirname(__file__), "data")


def test_real_abif_file():
    data = open(os.path.join(D, "biopython_3730.ab1"), "rb").read()
    tags = read_abif(data)
    assert tags["FWO_1"].startswith("GATC") and len(tags["PBAS2"]) == 1165
    tr = load_trace(data, "3730")
    a, b = tr.trim
    assert 0 < a < 30 and b > 1000                       # low-quality ends trimmed
    assert sum(tr.mixed[a:b]) <= 2 and tr.het_indel_start is None   # clean homozygous read


def _run(tmp_path, files):
    ex = write_examples(str(tmp_path / "in"))
    ann = json.load(open(ex["annotation"]))
    res = Pipeline(files, [ex["reference_dna"]], annotation=ann).run()
    return {s["query"].id: s for s in res["samples"]}


def test_het_snv_and_het_deletion_from_trace(tmp_path):
    a2 = HBB_CDS[:19] + "T" + HBB_CDS[20:123] + HBB_CDS[127:]
    f = str(tmp_path / "p1.ab1")
    write_synthetic_ab1(HBB_CDS, a2, f, name="p1")
    s = _run(tmp_path, [f])["p1"]
    got = {v["hgvs_c"]: v for v in s["variants"]}
    assert set(got) == {"c.20A>T", "c.126_129del"}
    assert all(v["zygosity"] == "heterozygous" for v in got.values())
    assert all(v.get("trace_svg") for v in got.values())


def test_reverse_read_and_homozygous(tmp_path):
    f = str(tmp_path / "p2.ab1")
    alt = HBB_CDS[:117] + "T" + HBB_CDS[118:]
    write_synthetic_ab1(reverse_complement(alt), None, f, name="p2")          # homozygous, reverse read
    s = _run(tmp_path, [f])["p2"]
    v = s["variants"][0]
    assert v["hgvs_c"] == "c.118C>T" and v["zygosity"].startswith("homozygous")
    assert any("reverse" in n for n in s["qc"]["notes"])


def test_het_deletion_in_reverse_read(tmp_path):
    """V11: a heterozygous deletion read on the reverse strand. Turned round, the double peaks lie BEFORE
    the deletion; they must be resolved as one heterozygous deletion, not called as many SNVs."""
    a2 = HBB_CDS[:123] + HBB_CDS[127:]                                       # c.126_129del on one allele
    f = str(tmp_path / "p3.ab1")
    write_synthetic_ab1(reverse_complement(HBB_CDS), reverse_complement(a2), f, name="p3")
    s = _run(tmp_path, [f])["p3"]
    got = {v["hgvs_c"]: v for v in s["variants"]}
    assert set(got) == {"c.126_129del"}, sorted(got)
    assert got["c.126_129del"]["zygosity"] == "heterozygous"
