import json
import os

import pytest

from clavagen.cohort import read_patient_table, run_cohort
from clavagen.examples import HBB_CDS, write_examples


def test_cohort_run(tmp_path):
    ex = write_examples(str(tmp_path / "in"))
    pts = [
        {"id": "P-001", "sex": "female", "age": "24", "phenotypes": "microcytic anaemia",
         "evidence": "c.20A>T=PS3,PS4", "query": HBB_CDS[:19] + "W" + HBB_CDS[20:]},
        {"id": "P-002", "sex": "male", "query": ">p2\n" + HBB_CDS[:117] + "T" + HBB_CDS[118:]},
        {"id": "P-003", "query": HBB_CDS},                       # no variants
    ]
    ann = json.load(open(ex["annotation"]))
    summ = run_cohort(pts, ex["reference_dna"], str(tmp_path / "out"), annotation=ann, online=False)
    folder = tmp_path / "out" / summ["name"]
    assert [e["ok"] for e in summ["patients"]] == [True, True, True]
    for pid in ("P-001", "P-002", "P-003"):
        assert (folder / pid / "report.html").exists()
    tsv = (folder / "cohort_summary.tsv").read_text()
    assert "P-001" in tsv and "c.20A>T" in tsv and "P-002" in tsv and "No variants detected" in tsv
    p1 = json.load(open(folder / "P-001" / "clinical_summary.json"))
    v = p1["samples"][0]["variants"][0]
    assert "PS3" in v["criteria_met"] and p1["meta"]["patient"]["sex"] == "female"
    p2 = json.load(open(folder / "P-002" / "clinical_summary.json"))
    assert "PS3" not in p2["samples"][0]["variants"][0]["criteria_met"]   # evidence stays per patient
    assert (folder / "cohort_summary.html").exists()


def test_duplicate_ids_and_missing_sequence(tmp_path):
    with pytest.raises(ValueError):
        run_cohort([{"id": "A", "query": "ACGT"}, {"id": "A", "query": "ACGT"}], HBB_CDS, str(tmp_path))
    with pytest.raises(ValueError):
        run_cohort([{"id": "A", "query": ""}], HBB_CDS, str(tmp_path))


def test_patient_table_csv(tmp_path):
    fa = tmp_path / "p2.fa"
    fa.write_text(">p2\n" + HBB_CDS)
    csvp = tmp_path / "patients.csv"
    csvp.write_text('id,sex,age,indication,phenotypes,evidence,query\n'
                    f'P1,female,30,Screening,"anaemia, splenomegaly",c.20A>T=PS3,{HBB_CDS}\n'
                    'P2,male,5,,,,p2.fa\n')
    pts = read_patient_table(str(csvp))
    assert pts[0]["phenotypes"] == ["anaemia", "splenomegaly"] and pts[0]["evidence"] == ["c.20A>T=PS3"]
    assert pts[1]["query"].startswith(">p2")
