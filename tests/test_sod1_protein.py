"""Regression test for the user-reported SOD1 case: protein input, NP_ reference header,
two ClinVar 2-star pathogenic variants (p.Lys4Glu, p.Ala5Ser). Online services are simulated."""

import os

import pytest

from clavagen import online
from clavagen.pipeline import Pipeline

D = os.path.join(os.path.dirname(__file__), "data")
CDS = ("ATGGCGACGAAGGCCGTGTGCGTGCTGAAGGGCGACGGCCCAGTGCAGGGCATCATCAATTTCGAGCAGAAGGAAAGTAATGGACCAGTGAAGGTGTGGGGAAGC"
       "ATTAAAGGACTGACTGAAGGCCTGCATGGATTCCATGTTCATGAGTTTGGAGATAATACAGCAGGCTGTACCAGTGCAGGTCCTCACTTTAATCCTCTATCCAG"
       "AAAACACGGTGGGCCAAAGGATGAAGAGAGGCATGTTGGAGACTTGGGCAATGTGACTGCTGACAAAGATGGTGTGGCCGATGTGTCTATTGAAGATTCTGTG"
       "ATCTCACTCTCAGGAGACCATTGCATCATTGGCCGCACACTGGTGGTCCATGAAAAAGCAGATGACTTGGGCAAAGGTGGAAATGAAGAAAGTACAAAGACAG"
       "GAAACGCTGGAAGTCGTTTGGCTTGTGGTGTAATTGGGATCGCCCAATAA")
UTR5, UTR3 = "GTTTGGGGCCAGAGTGGGCGAGGCGCGGAGGTCTGGCCTATAAAGTAGTCGCGGAGACGGGGTGCTGGTTTGCGTCGTAGTCTCCTGCAGCGTCTGGGGTTTCCGTTGCAGTCCTCGGAACCAGGACCTCGGCGTGGCCTAGCGAGTT", "GCAGATGACTTGGGCAAAGGTGGAAATGAAGAAAGTACAAAGACAGGAAAC"
MRNA = UTR5 + CDS + UTR3
CS = len(UTR5) + 1
CE = len(UTR5) + len(CDS)
MULTI = "criteria provided, multiple submitters, no conflicts"
CLINVAR = {
    "101": ("NM_000454.5(SOD1):c.10A>G (p.Lys4Glu)", "Pathogenic", MULTI),
    "102": ("NM_000454.5(SOD1):c.13G>T (p.Ala5Ser)", "Pathogenic", MULTI),
    "103": ("NM_000454.5(SOD1):c.14C>T (p.Ala5Val)", "Pathogenic", "reviewed by expert panel"),
    "104": ("NM_000454.5(SOD1):c.14C>A (p.Ala5Asp)", "Likely pathogenic", "criteria provided, single submitter"),
    "105": ("NM_000454.5(SOD1):c.272A>C (p.Asp91Ala)", "Pathogenic", MULTI),
    "106": ("NM_000454.5(SOD1):c.20G>C (p.Cys7Ser)", "Pathogenic", "criteria provided, single submitter"),
}


def _gb():
    seq = MRNA.lower()
    origin = "\n".join(f"{i + 1:>9} " + " ".join(seq[j:j + 10] for j in range(i, min(i + 60, len(seq)), 10))
                       for i in range(0, len(seq), 60))
    L = len(MRNA)
    e1 = CS - 1 + 72
    return (f"LOCUS       NM_000454    {L} bp    mRNA\nVERSION     NM_000454.5\nFEATURES             Location/Qualifiers\n"
            f"     gene            1..{L}\n                     /gene=\"SOD1\"\n     exon            1..{e1}\n"
            f"     CDS             {CS}..{CE}\n                     /protein_id=\"NP_000445.1\"\n"
            f"     exon            {e1 + 1}..{e1 + 97}\n     exon            {e1 + 98}..{e1 + 167}\n"
            f"     exon            {e1 + 168}..{L}\nORIGIN\n{origin}\n//\n")


def fake_http(url, data=None, headers=None, timeout=30, use_cache=True, kind="json"):
    if "efetch.fcgi" in url and "db=protein" in url:
        return ('LOCUS       NP_000445    154 aa\nVERSION     NP_000445.1\nFEATURES\n     gene            1..154\n'
                '                     /gene="SOD1"\n     CDS             1..154\n'
                '                     /coded_by="NM_000454.5:156..620"\nORIGIN\n//\n')
    if "efetch.fcgi" in url:
        return _gb()
    if "esearch.fcgi" in url:
        assert "SOD1" in url
        return {"esearchresult": {"idlist": list(CLINVAR)}}
    if "esummary.fcgi" in url:
        res = {"uids": list(CLINVAR)}
        for k, (t, c, r) in CLINVAR.items():
            res[k] = {"title": t, "germline_classification": {
                "description": c, "review_status": r,
                "trait_set": [{"trait_name": "Amyotrophic lateral sclerosis type 1"}]}}
        return {"result": res}
    if "uniprot" in url:
        return {"results": [{"primaryAccession": "P00441", "features": [], "comments": [
            {"commentType": "DISEASE", "disease": {"diseaseId": "Amyotrophic lateral sclerosis 1",
             "description": "A neurodegenerative disorder ... inherited in an autosomal dominant or, rarely, "
                            "autosomal recessive manner.", "diseaseCrossReference": {"database": "MIM", "id": "105400"}}}]}]}
    if "gnomad" in url:
        return {"data": {"gene": {"gnomad_constraint": {"oe_lof_upper": 1.1, "mis_z": 1.5}}}}
    if "rest.ensembl.org/vep" in url:
        out = []
        for h in data["hgvs_notations"]:
            rec = {"input": h, "transcript_consequences": [{"transcript_id": "NM_000454.5", "revel": 0.9}]}
            out.append(rec)   # absent from gnomAD
        return out
    raise AssertionError(url)


@pytest.fixture
def mocked(monkeypatch):
    monkeypatch.setattr(online, "http", fake_http)
    online._dead_hosts.clear()


def _run(**kw):
    res = Pipeline([os.path.join(D, "patient1_SOD1.fa"), os.path.join(D, "patient2_SOD1.fa")],
                   [os.path.join(D, "Ref_SOD1.fa")], online=True, pairs=None, **kw).run()
    return {s["query"].id: s["variants"][0] for s in res["samples"]}, res


def test_np_header_resolves_gene_and_classifies(mocked):
    vs, res = _run()
    k4e, a5s = vs["Patient1"], vs["Patient2"]
    assert k4e["hgvs_p"] == "p.(Lys4Glu)" and a5s["hgvs_p"] == "p.(Ala5Ser)"
    assert res["samples"][0]["context"]["gene"] == "SOD1"
    # nucleotide change inferred from the RefSeq codon
    assert k4e["hgvs_c"] == "c.10A>G" and k4e["nucleotide_inferred"]
    assert a5s["hgvs_c"] == "c.13G>T"
    # reported classification comes from the 2-star ClinVar record
    assert k4e["classification"] == "Pathogenic" and k4e["clinvar_stars"] == 2
    assert a5s["classification"] == "Pathogenic"
    assert "ClinVar" in k4e["classification_basis"]
    # automated ACMG evidence now available: PM2 (absent from gnomAD), PP3 (REVEL), PM5/PM1 at residue 5
    assert "PM2_Supporting" in k4e["criteria_met"] and any(c.startswith("PP3") for c in k4e["criteria_met"])
    assert "PM5" in a5s["criteria_met"] or "PM1" in a5s["criteria_met"]
    assert res["samples"][0]["context"]["inheritance"] in ("AD", "AD/AR", "AD (rarely AR)")


def test_acmg_only_mode(mocked):
    vs, _ = _run(use_clinvar=False)
    assert vs["Patient1"]["classification"] == vs["Patient1"]["acmg_classification"]
    assert vs["Patient1"]["clinvar_classification"] == "Pathogenic"


def test_ssl_error_message():
    import ssl
    import urllib.error
    msg = online.explain_error(urllib.error.URLError(ssl.SSLCertVerificationError("CERTIFICATE_VERIFY_FAILED")))
    assert "Install Certificates" in msg


def _run4(**kw):
    res = Pipeline([os.path.join(D, "patient4_SOD1.fa")], [os.path.join(D, "Ref_SOD1.fa")], online=True, **kw).run()
    return res["samples"][0]["variants"][0]


def test_one_star_clinvar_is_used_by_default(mocked):
    v = _run4()
    assert v["hgvs_p"] == "p.(Cys7Ser)"
    assert sorted(v["hgvs_c_candidates"]) == ["c.19T>A", "c.20G>C"]       # TGC -> AGC or TCC
    assert v["clinvar_stars"] == 1 and v["classification"] == "Pathogenic"
    assert "limited review" in v["classification_basis"]


def test_min_stars_setting(mocked):
    v = _run4(clinvar_min_stars=2)
    assert v["classification"] == v["acmg_classification"]              # 1-star record not used
    assert v["clinvar_classification"] == "Pathogenic"                   # but still shown
