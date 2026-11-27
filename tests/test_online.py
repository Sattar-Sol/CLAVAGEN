"""Online annotation tests with simulated service responses (no network needed)."""

import json

import pytest

from clavagen import online
from clavagen.examples import HBB_CDS
from clavagen.pipeline import Pipeline

UTR5 = "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACC"  # 50 nt
UTR3 = "GCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGG"
MRNA = UTR5 + HBB_CDS + UTR3


def _genbank():
    seq = MRNA.lower()
    origin = "\n".join(f"{i + 1:>9} " + " ".join(seq[j:j + 10] for j in range(i, min(i + 60, len(seq)), 10))
                       for i in range(0, len(seq), 60))
    L = len(MRNA)
    return f"""LOCUS       NM_000518    {L} bp    mRNA    linear   PRI
DEFINITION  Homo sapiens hemoglobin subunit beta (HBB), mRNA.
VERSION     NM_000518.5
FEATURES             Location/Qualifiers
     source          1..{L}
     gene            1..{L}
                     /gene="HBB"
     exon            1..142
                     /gene="HBB"
     CDS             51..494
                     /gene="HBB"
                     /protein_id="NP_000509.1"
     exon            143..365
                     /gene="HBB"
     exon            366..{L}
                     /gene="HBB"
ORIGIN
{origin}
//
"""


CLINVAR = {
    "1": ("NM_000518.5(HBB):c.118C>T (p.Gln40Ter)", "Pathogenic"),
    "2": ("NM_000518.5(HBB):c.20A>T (p.Glu7Val)", "Pathogenic"),
    "3": ("NM_000518.5(HBB):c.19G>A (p.Glu7Lys)", "Pathogenic"),
    "4": ("NM_000518.5(HBB):c.126_129del (p.Phe42fs)", "Pathogenic"),
    "5": ("NM_000518.5(HBB):c.92+1G>A", "Pathogenic/Likely pathogenic"),
    "6": ("NM_000518.5(HBB):c.9T>C (p.His3=)", "Benign"),
    "7": ("NM_000518.5(HBB):c.20A>G (p.Glu7Gly)", "Likely pathogenic"),
}


def fake_http(url, data=None, headers=None, timeout=30, use_cache=True, kind="json"):
    if "esearch.fcgi" in url:
        return {"esearchresult": {"idlist": list(CLINVAR)}}
    if "esummary.fcgi" in url:
        res = {"uids": list(CLINVAR)}
        for k, (t, c) in CLINVAR.items():
            res[k] = {"uid": k, "title": t, "germline_classification": {"description": c,
                                                                       "review_status": "criteria provided"}}
        return {"result": res}
    if "efetch.fcgi" in url:
        return _genbank()
    if "uniprot" in url:
        return {"results": [{"primaryAccession": "P68871", "features": [
            {"type": "Domain", "location": {"start": {"value": 3}, "end": {"value": 147}}, "description": "Globin"},
            {"type": "Binding site", "location": {"start": {"value": 64}, "end": {"value": 64}}},
            {"type": "Binding site", "location": {"start": {"value": 93}, "end": {"value": 93}}}],
            "comments": [{"commentType": "FUNCTION", "texts": [{"value": "Involved in oxygen transport."}]},
                         {"commentType": "DISEASE", "disease": {
                             "diseaseId": "Sickle cell disease",
                             "description": "An autosomal recessive haemoglobinopathy...",
                             "diseaseCrossReference": {"database": "MIM", "id": "603903"}}}]}]}
    if "gnomad" in url:
        return {"data": {"gene": {"symbol": "HBB", "gnomad_constraint": {"oe_lof_upper": 1.2, "mis_z": 0.5}}}}
    if "rest.ensembl.org/vep" in url:
        out = []
        for h in data["hgvs_notations"]:
            rec = {"input": h, "transcript_consequences": [{"transcript_id": "NM_000518.5"}]}
            if h.endswith("c.9T>C"):
                rec["colocated_variants"] = [{"id": "rs713040", "frequencies": {"C": {"gnomade": 0.21, "gnomadg": 0.25}}}]
            if h.endswith("c.20A>T"):
                rec["colocated_variants"] = [{"id": "rs334", "clin_sig": ["pathogenic"],
                                              "frequencies": {"T": {"gnomade": 0.004}}}]
                rec["transcript_consequences"][0].update({"revel": 0.71, "spliceai": {"DS_AG": 0.0, "DS_DL": 0.02}})
            out.append(rec)
        return out
    raise AssertionError(url)


@pytest.fixture
def mocked(monkeypatch):
    monkeypatch.setattr(online, "http", fake_http)


def test_parsers(mocked):
    gb = online.fetch_refseq("NM_000518.5")
    assert gb["cds"] == [51, 494] and gb["exon_boundaries_c"] == [92, 315] and gb["gene"] == "HBB"
    kv = online.fetch_clinvar("HBB")
    assert {k["hgvs_c"] for k in kv} == {"c.118C>T", "c.20A>T", "c.19G>A", "c.126_129del", "c.92+1G>A", "c.20A>G"}
    ann = online.auto_annotation("HBB")
    assert ann["transcript"] == "NM_000518.5" and ann["lof_mechanism"] is True
    assert ann["inheritance"] == "AR" and 7 in ann["hotspots"]


def test_manual_known_parser():
    t = "c.20A>T | p.Glu7Val | Pathogenic | Lab database 2024\nNM_000518.5:c.118C>T, LP\n# comment\nE7K, pathogenic"
    kv = online.parse_known_variants_text(t)
    assert kv[0] == {"hgvs_c": "c.20A>T", "hgvs_p": "p.Glu7Val", "classification": "Pathogenic", "source": "Lab database 2024"}
    assert kv[1]["hgvs_c"] == "c.118C>T" and kv[1]["classification"] == "Likely pathogenic"
    assert kv[2]["hgvs_p"] == "p.Glu7Lys"


def test_pipeline_online_with_accession_reference(mocked):
    q = HBB_CDS[:8] + "C" + HBB_CDS[9:117] + "T" + HBB_CDS[118:]   # c.9T>C + c.118C>T
    pl = Pipeline([f"P2={q}"], ["NM_000518.5"], online=True,
                  known_variants=[{"hgvs_c": "c.300G>A", "hgvs_p": "p.Trp100Ter", "classification": "Pathogenic",
                                   "source": "In-house"}])
    res = pl.run()
    s = res["samples"][0]
    by = {v["hgvs_c"]: v for v in s["variants"]}
    assert set(by) == {"c.9T>C", "c.118C>T"}
    assert by["c.9T>C"]["classification"] in ("Benign", "Likely benign")       # BA1 from gnomAD AF 0.25
    assert by["c.118C>T"]["classification"] in ("Pathogenic", "Likely pathogenic")
    assert any("ClinVar" in a for a in by["c.118C>T"]["database_assertions"])
    assert s["context"]["cds"] == [51, 494]
    assert any(k.get("source") == "In-house" for k in s["context"]["known_variants"])


def test_offline_failure_is_graceful(monkeypatch):
    def boom(*a, **k):
        raise online.OnlineError("network unreachable")
    monkeypatch.setattr(online, "http", boom)
    res = Pipeline([f"P={HBB_CDS[:117]}T{HBB_CDS[118:]}"], [f"HBB_NM_000518.5={HBB_CDS}"], gene="HBB", online=True).run()
    notes = " ".join(res["samples"][0]["qc"]["notes"])
    assert "network unreachable" in notes
    assert res["samples"][0]["variants"][0]["consequence"] == "stop_gained"


# ---------------------------------------------------------------- 1.0.2: old transcript versions & strand
def _gb(acc, seq, cds, exons):
    return {"transcript": acc, "seq": seq, "cds": list(cds), "exon_boundaries_c": exons, "gene": "TEST"}


def test_old_transcript_version_mapped_to_current(monkeypatch):
    from clavagen import online
    seq = "GGG" + "ATGAAACCCTTTGGGTAA" + "TTT"
    recs = {"NM_000001": _gb("NM_000001.4", "AAAA" + seq, (8, 25), [6]),     # longer 5'UTR, same CDS
            "NM_000001.3": _gb("NM_000001.3", seq, (4, 21), [6])}
    monkeypatch.setattr(online, "fetch_refseq", lambda a: recs[a])
    acc, note = online.current_transcript_version("NM_000001.3")
    assert acc == "NM_000001.4" and "not the current RefSeq version" in note
    # same version -> unchanged, no note
    recs["NM_000001"] = recs["NM_000001.3"]
    assert online.current_transcript_version("NM_000001.3") == ("NM_000001.3", None)


def test_old_transcript_reference_base_change(monkeypatch):
    """NM_000277.1 vs .3 (PAH): same CDS length, two reference bases differ, no exon annotation in .1."""
    from clavagen import online
    recs = {"NM_000002": _gb("NM_000002.2", "ATGAAACCCTAA", (1, 12), [6]),
            "NM_000002.1": _gb("NM_000002.1", "ATGAAAGCCTAA", (1, 12), [])}
    monkeypatch.setattr(online, "fetch_refseq", lambda a: recs[a])
    m = online.transcript_mapping("NM_000002.1")
    assert m["acc"] == "NM_000002.2" and m["diff"] == {7} and not m["exons_ok"]
    assert online.mappable("c.4A>G", m) and online.mappable("c.2_3del", m)
    assert not online.mappable("c.7G>T", m)          # reference base changed there
    assert not online.mappable("c.6+1G>A", m)        # exon boundaries not comparable
    assert not online.mappable("c.-5A>G", m)         # UTR numbering can shift


def test_old_transcript_with_different_cds_length_not_mapped(monkeypatch):
    from clavagen import online
    recs = {"NM_000003": _gb("NM_000003.2", "ATGAAACCCGGGTAA", (1, 15), []),
            "NM_000003.1": _gb("NM_000003.1", "ATGAAACCCTAA", (1, 12), [])}
    monkeypatch.setattr(online, "fetch_refseq", lambda a: recs[a])
    acc, note = online.current_transcript_version("NM_000003.1")
    assert acc == "NM_000003.1" and "differs" in note


def test_revel_added_from_ensembl_transcript(monkeypatch):
    """VEP gives REVEL only on Ensembl transcripts: an HGVS lookup on RefSeq never has it."""
    from clavagen import online
    calls = []

    def fake(url, data=None, **k):
        if "/hgvs" in url:
            return [{"input": h, "seq_region_name": "12", "start": 102843691, "allele_string": "T/C", "strand": -1,
                     "transcript_consequences": [{"transcript_id": "NM_000277.3"}]} for h in data["hgvs_notations"]]
        calls.append(data["variants"])
        return [{"input": "12 102843691 . A G . . .", "transcript_consequences": [
            {"transcript_id": "ENST00000307000", "revel": 0.5},
            {"transcript_id": "ENST00000553106", "mane_select": "NM_000277.3", "revel": 0.817}]}]
    monkeypatch.setattr(online, "http", fake)
    res = online.vep_lookup(["NM_000277.3:c.1154T>C"], "NM_000277.3")
    assert calls == [["12 102843691 . A G . . ."]]         # minus strand: T/C on the transcript = A/G genomic
    assert res["NM_000277.3:c.1154T>C"]["REVEL"] == 0.817 and "_loc" not in res["NM_000277.3:c.1154T>C"]


def test_vep_alt_allele_from_allele_string():
    from clavagen.online import _vep_alt, parse_vep_record
    rec = {"allele_string": "G/A", "colocated_variants": [
        {"frequencies": {"A": {"gnomade": 0.03}, "T": {"gnomade": 0.9}}}]}
    assert _vep_alt(rec) == "A"
    assert parse_vep_record(rec, _vep_alt(rec))["af"] == 0.03
    assert _vep_alt({"allele_string": "GA/-"}) is None


def test_vep_batch_survives_one_bad_hgvs(monkeypatch):
    from clavagen import online
    def fake(chunk):
        if any("BAD" in h for h in chunk):
            raise online.OnlineError("rest.ensembl.org: HTTP 400")
        return [{"input": h, "allele_string": "C/T",
                 "colocated_variants": [{"frequencies": {"T": {"gnomadg": 0.2}}}]} for h in chunk]
    monkeypatch.setattr(online, "_vep_post", fake)
    res = online.vep_lookup(["NM_1.1:c.1C>T", "NM_1.1:c.BAD", "NM_1.1:c.3C>T"])
    assert set(res) == {"NM_1.1:c.1C>T", "NM_1.1:c.3C>T"} and res["NM_1.1:c.1C>T"]["af"] == 0.2


def test_delins_with_deleted_bases_matches():
    from clavagen.annotation import normalize_hgvs
    assert normalize_hgvs("c.470_471delGAinsAC") == normalize_hgvs("c.470_471delinsAC")
