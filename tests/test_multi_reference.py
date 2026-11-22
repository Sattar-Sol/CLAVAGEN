import json

from clavagen.annotation import user_context
from clavagen.cohort import run_cohort
from clavagen.examples import HBB_CDS, SOD1_CDS

IU = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M"}


def var(seq, pos, alt, het):
    s = list(seq)
    s[pos - 1] = IU["".join(sorted(s[pos - 1] + alt))] if het else alt
    return "".join(s)


REFS = [{"fasta": f">NM_000518.5 (HBB) beta globin\n{HBB_CDS}", "id": "NM_000518.5"},
        {"fasta": f">NM_000454.5 (SOD1) superoxide dismutase\n{SOD1_CDS}", "id": "NM_000454.5"}]
ANN = {"references": {
    "NM_000518.5": {"gene": "HBB", "known_variants": [{"hgvs_c": "c.20A>T", "hgvs_p": "p.Glu7Val",
                                                       "classification": "Pathogenic", "source": "test list"}]},
    "NM_000454.5": {"gene": "SOD1", "inheritance": "AD/AR",
                    "known_variants": [{"hgvs_c": "c.272A>C", "hgvs_p": "p.Asp91Ala",
                                        "classification": "Pathogenic", "source": "test list"}]}}}


def test_two_genes_one_run_auto_matching(tmp_path):
    pts = [{"id": "H1", "query": var(HBB_CDS, 20, "T", True)},
           {"id": "S1", "query": var(SOD1_CDS, 272, "C", False)},
           {"id": "S2", "query": var(SOD1_CDS, 272, "C", True)}]
    summ = run_cohort(pts, REFS, str(tmp_path), annotation=ANN, online=False)
    out = {}
    for e in summ["patients"]:
        js = json.load(open(tmp_path / summ["name"] / e["folder"] / "clinical_summary.json"))
        smp = js["samples"][0]
        out[e["id"]] = (smp["reference"]["id"], smp["variants"][0]["gene"], smp["variants"][0]["classification"],
                        smp["clinical_summary"]["status"])
    assert out["H1"][:3] == ("NM_000518.5", "HBB", "Pathogenic")
    assert out["S1"][:3] == ("NM_000454.5", "SOD1", "Pathogenic") and "biallelic" in out["S1"][3].lower()
    assert "carrier" in out["S2"][3].lower()          # per-reference AD/AR setting applied to SOD1 only
    assert "carrier" in out["H1"][3].lower()


def test_patient_pinned_to_reference(tmp_path):
    pts = [{"id": "H1", "reference": "HBB", "query": var(HBB_CDS, 20, "T", True)},
           {"id": "X", "reference": "NM_999999.1", "query": HBB_CDS}]
    summ = run_cohort(pts, REFS, str(tmp_path), annotation=ANN, online=False)
    ok = {e["id"]: e["ok"] for e in summ["patients"]}
    assert ok == {"H1": True, "X": False}
    assert "not loaded" in next(e["error"] for e in summ["patients"] if e["id"] == "X")


def test_settings_match_without_version():
    assert user_context({"references": {"NM_000454": {"inheritance": "AD/AR"}}}, "NM_000454.5")["inheritance"] == "AD/AR"


def test_window_payload_with_reference_rows(tmp_path):
    from clavagen.webui import parse_reference_text, run_payload
    rows = parse_reference_text(REFS[0]["fasta"] + "\n" + REFS[1]["fasta"])
    assert [(r["id"], r["gene"]) for r in rows] == [("NM_000518.5", "HBB"), ("NM_000454.5", "SOD1")]
    rows[1]["inheritance"] = "AD/AR"
    acc = parse_reference_text("NM_000518.5\nNM_000454")
    assert acc[0]["accession"] == "NM_000518.5" and acc[1]["id"] == "NM_000454"
    d = run_payload({"references": rows, "online": "0", "patients": [{"id": "S2", "query": var(SOD1_CDS, 272, "C", True)}],
                     "known": "c.272A>C | p.Asp91Ala | Pathogenic | lab"}, str(tmp_path))
    p = d["patients"][0]
    assert p["samples"][0]["variants"][0]["cls"] == "Pathogenic"
    assert "carrier" in p["samples"][0]["status"].lower()
