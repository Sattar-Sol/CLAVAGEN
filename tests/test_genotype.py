"""V12: patient-level genotype across files, phase from the family, variant types, VUS included."""
import json

from clavagen.cohort import run_cohort
from clavagen.examples import HBB_CDS, write_examples
from clavagen.genetic_code import translate

IU = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M"}


def het(seq, *changes):
    s = list(seq)
    for pos, alt in changes:
        s[pos - 1] = IU["".join(sorted(s[pos - 1] + alt))]
    return "".join(s)


def amplicons(seq):
    """Split a sequence into two overlapping 'amplicon' files, as a lab would."""
    return [f">ampA\n{seq[:250]}", f">ampB\n{seq[200:]}"]


def _run(tmp_path, pts):
    ex = write_examples(str(tmp_path / "in"))
    ann = json.load(open(ex["annotation"]))
    summ = run_cohort(pts, ex["reference_dna"], str(tmp_path / "out"), annotation=ann, online=False)
    return tmp_path / "out" / summ["name"]


def _clin(folder, pid):
    js = json.load(open(folder / pid / "clinical_summary.json"))
    return [s["clinical_summary"] for s in js["samples"]]


# c.118C>T p.(Gln40Ter) is P in the demo annotation; c.250G>A and c.300T>A are unannotated (VUS)
VUS_A, VUS_B = (250, "A"), (300, "A")


def test_compound_het_vus_across_two_files(tmp_path):
    pts = [
        {"id": "Kid", "family": "F", "relationship": "proband", "affected": "affected",
         "query": amplicons(het(HBB_CDS, VUS_A, VUS_B))},
        {"id": "Mum", "family": "F", "relationship": "mother", "affected": "unaffected",
         "query": het(HBB_CDS, VUS_A)},
        {"id": "Dad", "family": "F", "relationship": "father", "affected": "unaffected",
         "query": het(HBB_CDS, VUS_B)},
    ]
    folder = _run(tmp_path, pts)
    for c in _clin(folder, "Kid"):                       # every file of the patient gives the same genotype
        assert c["genotype"]["genotype"] == "compound heterozygous" and c["genotype"]["phase"] == "trans"
        assert c["status"].startswith("Compound heterozygous, in trans (confirmed")
        assert "both variants VUS" in c["status"]
        assert "missense / missense" in c["genotype"]["allele_types"]


def test_compound_het_null_missense_phase_unknown(tmp_path):
    pts = [{"id": "Solo", "query": amplicons(het(HBB_CDS, (118, "T"), VUS_A))}]   # nonsense + VUS, no parents
    folder = _run(tmp_path, pts)
    c = _clin(folder, "Solo")[0]
    g = c["genotype"]
    assert g["phase"] == "unknown" and "possible compound heterozygote" in g["genotype"]
    assert "null (nonsense)" in g["allele_types"] and "missense" in g["allele_types"]
    assert any("phase" in r.lower() for r in c["recommendations"])


def test_protein_alleles_on_separate_sequences(tmp_path):
    p = translate(HBB_CDS, True)
    a1 = p[:6] + "V" + p[7:]                              # p.Glu7Val
    a2 = p[:39] + "*"                                     # p.Gln40Ter
    folder = _run(tmp_path, [{"id": "Prot", "query": f">allele1\n{a1}\n>allele2\n{a2}"}])
    c = _clin(folder, "Prot")[0]
    assert "zygosity not determinable" in c["genotype"]["genotype"]
    assert any("separate input sequences" in x for x in c["interpretation"])
    assert "possible compound heterozygote" in c["status"]


def test_consensus_indel_not_called_homozygous(tmp_path):
    """A heterozygous indel cannot be shown in a text/consensus sequence: it must not be reported as a
    homozygous (biallelic) genotype - zygosity is 'unresolved' until a trace is used."""
    seq = HBB_CDS[:125] + HBB_CDS[129:]                       # c.126_129del written into a consensus
    folder = _run(tmp_path, [{"id": "Cons", "query": seq}])
    c = _clin(folder, "Cons")[0]
    assert c["genotype"]["genotype"] == "indel, zygosity unresolved"
    assert "biallelic" not in c["status"].lower() and "zygosity unresolved" in c["status"]
    assert any(".ab1" in r for r in c["recommendations"])


def test_offline_gene_from_header_and_gene_keyed_setting(tmp_path):
    """Offline, without an annotation file: the gene symbol is read from the FASTA header ('gene=...')
    and a per-reference setting keyed by that gene symbol is applied."""
    from acmg_varcall.cli import _apply_ref_settings  # old name keeps working too
    ref = tmp_path / "ref.fa"
    ref.write_text(">NM_TEST.1 gene=MYGENE test transcript\n" + HBB_CDS + "\n")
    ann = _apply_ref_settings({}, ["MYGENE:inheritance=AR"])
    pts = [{"id": "Kid", "family": "F", "relationship": "proband", "affected": "affected",
            "query": amplicons(het(HBB_CDS, VUS_A, VUS_B))},
           {"id": "Mum", "family": "F", "relationship": "mother", "query": het(HBB_CDS, VUS_A)},
           {"id": "Dad", "family": "F", "relationship": "father", "query": het(HBB_CDS, VUS_B)}]
    summ = run_cohort(pts, str(ref), str(tmp_path / "out"), annotation=ann, online=False)
    c = _clin(tmp_path / "out" / summ["name"], "Kid")[0]
    assert c["gene"] == "MYGENE" and c["inheritance"] == "AR"
    assert c["genotype"]["genotype"] == "compound heterozygous"
