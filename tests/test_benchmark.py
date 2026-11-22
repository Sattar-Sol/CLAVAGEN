import json
import os

import pytest

from clavagen import online
from clavagen.benchmark import load_erepo, run_benchmark

from test_online import fake_http   # simulated NCBI / UniProt / gnomAD / VEP for HBB

EREPO = "\t".join(["#Variation", "ClinVar Variation Id", "Allele Registry Id", "HGVS Expressions", "HGNC Gene Symbol",
                   "Disease", "Mondo Id", "Mode of Inheritance", "Assertion", "Applied Evidence Codes (Met)",
                   "Applied Evidence Codes (Not Met)", "Expert Panel", "Retracted"]) + "\n"
ROWS = [
    ("NM_000518.5(HBB):c.118C>T (p.Gln40Ter)", "Pathogenic", "PVS1, PM3_VeryStrong, PM2_Supporting", "false"),
    ("NM_000518.5(HBB):c.20A>T (p.Glu7Val)", "Pathogenic", "PS3, PS4, PM2_Supporting", "false"),
    ("NM_000518.5(HBB):c.9T>C (p.His3=)", "Benign", "BA1", "false"),
    ("NM_000518.5(HBB):c.92+1G>A", "Pathogenic", "PVS1", "false"),               # intronic -> skipped
    ("NM_000518.5(HBB):c.19G>A (p.Glu7Lys)", "Pathogenic", "PS3", "true"),       # retracted -> skipped
]
for t, a, m, retr in ROWS:
    EREPO += "\t".join([t, "1", "CA1", t.split(" ")[0], "HBB", "Hb disorder", "MONDO:1", "AR", a, m, "", "Hemoglobinopathy VCEP",
                        retr]) + "\n"


@pytest.fixture
def mocked(monkeypatch):
    monkeypatch.setattr(online, "http", fake_http)
    online._dead_hosts.clear()


def test_load_erepo():
    recs, stats = load_erepo(EREPO)
    assert stats["used"] == 3 and stats["retracted"] == 1 and stats["intronic_or_complex"] == 1
    assert recs[0]["expert_codes"] == ["PM2", "PM3", "PVS1"]


def test_run_benchmark(mocked, tmp_path):
    recs, _ = load_erepo(EREPO)
    s = run_benchmark(recs, str(tmp_path), log=lambda *_: None)
    assert s["analysed"] == 3 and s["errors"] == 0
    rows = {l.split("\t")[2]: l.split("\t") for l in open(tmp_path / "benchmark_results.tsv").read().splitlines()[1:]}
    assert rows["c.9T>C"][4] in ("Benign", "Likely benign")
    assert rows["c.118C>T"][4] in ("Pathogenic", "Likely pathogenic")
    assert s["critical_discordance"] == 0
    assert "PVS1" in s["criteria"] and s["criteria"]["PVS1"]["tp"] == 1
    assert (tmp_path / "benchmark_report.html").exists()
