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
    import csv
    rows = {r["hgvs_c"]: r for r in csv.DictReader(open(tmp_path / "benchmark_results.tsv"), delimiter="\t")}
    assert rows["c.9T>C"]["tool_class"] in ("Benign", "Likely benign")
    assert rows["c.118C>T"]["tool_class"] in ("Pathogenic", "Likely pathogenic")
    assert rows["c.118C>T"]["variant_type"] == "nonsense" and rows["c.20A>T"]["variant_type"] == "missense"
    assert s["critical_discordance"] == 0
    assert "PVS1" in s["criteria"] and s["criteria"]["PVS1"]["tp"] == 1
    assert (tmp_path / "benchmark_report.html").exists()


# ---------------------------------------------------------------- CLAVAGEN 1.0.1: leakage and disagreement causes
def test_parse_criteria_keeps_strengths():
    from clavagen.benchmark import parse_criteria
    c = parse_criteria("PVS1, PM3_VeryStrong, PM2_Supporting; PS3_Moderate, BA1, PP3_Very Strong")
    assert c == {"PVS1": "VeryStrong", "PM3": "VeryStrong", "PM2": "Supporting", "PS3": "Moderate",
                 "BA1": "StandAlone", "PP3": "VeryStrong"}


def test_leave_one_out_recounts_hotspots_and_lof():
    """The tested variant must not count towards its own hotspot (PM1) or loss-of-function (PVS1) evidence."""
    from clavagen.online import derive_gene_evidence, recalc_gene_evidence
    kv = lambda c, p: {"hgvs_c": c, "hgvs_p": p, "classification": "Pathogenic"}  # noqa: E731
    known = [kv("c.10C>T", "p.(Arg4Cys)"), kv("c.11G>A", "p.(Arg4His)"),        # residue 4: exactly 2 -> hotspot
             kv("c.30C>A", "p.(Tyr10Ter)"), kv("c.40C>T", "p.(Gln14Ter)"), kv("c.50del", "p.(Lys17fs)")]  # 3 truncating
    ann = derive_gene_evidence({"known_variants": list(known)})
    assert ann["hotspots"] == [4] and ann["lof_mechanism"] is True
    # leave out the tested missense c.11G>A: residue 4 has one P/LP missense left -> no longer a hotspot
    ctx = dict(ann, known_variants=[k for k in known if k["hgvs_c"] != "c.11G>A"])
    recalc_gene_evidence(ctx)
    assert ctx["hotspots"] == [] and ctx["lof_mechanism"] is True
    # leave out the tested nonsense c.40C>T: two truncating variants left -> LOF mechanism not supported
    ctx = dict(ann, known_variants=[k for k in known if k["hgvs_c"] != "c.40C>T"])
    recalc_gene_evidence(ctx)
    assert not ctx["lof_mechanism"] and ctx["hotspots"] == [4]
    # ... unless LOF comes from somewhere else: gnomAD constraint, or the built-in gene list / user annotation
    ctx = dict(ann, constraint={"loeuf": 0.2}, known_variants=[k for k in known if k["hgvs_c"] != "c.40C>T"])
    recalc_gene_evidence(ctx)
    assert ctx["lof_mechanism"] is True
    ctx = dict(ann, known_variants=[k for k in known if k["hgvs_c"] != "c.40C>T"])
    recalc_gene_evidence(ctx, fixed_hotspots={4, 99}, fixed_lof=True)
    assert ctx["lof_mechanism"] is True and ctx["hotspots"] == [4, 99]


def test_disagreement_categories():
    from clavagen.benchmark import analyse_disagreement
    base = {"status": "ok", "expert_class": "Pathogenic"}
    # experts reach P only with functional and case data; from computable criteria alone they get LP, as the tool
    r = analyse_disagreement(dict(base, tool_class="Likely pathogenic",
                                  expert_criteria={"PVS1": "VeryStrong", "PM2": "Supporting", "PS3": "Strong"},
                                  tool_criteria={"PVS1": "VeryStrong", "PM2": "Supporting"}))
    assert r["disagreement"] == "missing evidence" and r["expert_computable_class"] == "Likely pathogenic"
    assert r["expert_noncomputable"] == ["PS3"]
    # the tool misses a computable criterion the experts applied
    r = analyse_disagreement(dict(base, tool_class="Uncertain significance",
                                  expert_criteria={"PVS1": "VeryStrong", "PM2": "Supporting", "PM1": "Moderate"},
                                  tool_criteria={"PM2": "Supporting", "PM1": "Moderate"}))
    assert r["disagreement"] == "criterion assignment" and r["missed_codes"] == ["PVS1"]
    # same computable criteria, different strength
    r = analyse_disagreement(dict(base, tool_class="Likely pathogenic",
                                  expert_criteria={"PVS1": "VeryStrong", "PM2": "Supporting", "PM5": "Moderate"},
                                  tool_criteria={"PVS1": "Strong", "PM2": "Supporting", "PM5": "Moderate"}))
    assert r["disagreement"] == "evidence strength" and r["strength_diffs"] == ["PVS1: tool Strong / expert VeryStrong"]
    # agreement
    r = analyse_disagreement(dict(base, tool_class="Pathogenic", expert_criteria={"PVS1": "VeryStrong"},
                                  tool_criteria={"PVS1": "VeryStrong", "PM2": "Supporting"}))
    assert r["disagreement"] == ""


def test_benchmark_report_has_disagreement_section(mocked, tmp_path):
    recs, _ = load_erepo(EREPO)
    s = run_benchmark(recs, str(tmp_path), log=lambda *_: None)
    assert set(s["disagreement"]) == {"missing evidence", "criterion assignment", "evidence strength"}
    assert sum(s["disagreement"].values()) == s["discordant"]
    assert s["computable_three_tier_agreement"] is not None
    html = open(tmp_path / "benchmark_report.html").read()
    assert "Why the classifications differ" in html and "recounted without it" in html
    head = open(tmp_path / "benchmark_results.tsv").readline().rstrip("\n").split("\t")
    assert "disagreement" in head and "expert_computable_class" in head


def test_leakage_fix_end_to_end(monkeypatch, tmp_path):
    """Full pipeline: a missense variant must not make its own residue a hotspot (PM1) in a leave-one-out run,
    while a genuinely different pathogenic missense at the same residue still gives PM5. The shared
    annotation cache must not carry the exclusion over to the next run."""
    from clavagen import pipeline as pl_mod
    from clavagen.examples import HBB_CDS
    from clavagen.online import derive_gene_evidence

    def fake_auto(gene, tx, log=None):
        kv = [{"hgvs_c": "c.199A>G", "hgvs_p": "p.(Lys67Glu)", "classification": "Pathogenic", "source": "ClinVar"},
              {"hgvs_c": "c.200A>C", "hgvs_p": "p.(Lys67Thr)", "classification": "Pathogenic", "source": "ClinVar"}]
        return derive_gene_evidence({"gene": "HBB", "known_variants": kv})

    monkeypatch.setattr(online, "auto_annotation", fake_auto)
    monkeypatch.setattr(online, "http", lambda *a, **k: (_ for _ in ()).throw(OSError("offline test")))
    ref = tmp_path / "ref.fa"
    ref.write_text(">HBB_test (HBB)\n" + HBB_CDS + "\n")
    query = HBB_CDS[:198] + "G" + HBB_CDS[199:]                       # c.199A>G p.(Lys67Glu)
    cache = {}

    def run(exclude):
        p = pl_mod.Pipeline([f"q={query}"], [str(ref)], gene="HBB", online=True, use_clinvar=False,
                            auto_cache=cache, log=lambda *_: None, exclude_known=exclude)
        v = p.run()["samples"][0]["variants"][0]
        met = {c["code"]: c["rationale"] for c in v["criteria"] if c["status"] == "met"}
        return v, met

    v, met = run({"c.199A>G"})                                          # leave-one-out
    assert v["hgvs_p"] == "p.(Lys67Glu)"
    assert not ("PM1" in met and "hotspot" in met["PM1"].lower()), met.get("PM1")
    assert "PM5" in met and "Lys67Thr" in met["PM5"]                   # the other variant is legitimate evidence
    v2, met2 = run(set())                                               # normal run, same cache
    assert "PM1" in met2 and "hotspot" in met2["PM1"].lower()           # hotspot back: cache was not altered


def test_protein_only_record_not_used_as_own_ps1(monkeypatch, tmp_path):
    """A ClinVar record of the tested variant that gives only the protein change must also be left out;
    otherwise it matches as 'same amino-acid change' and gives the variant PS1 from its own record."""
    from clavagen import pipeline as pl_mod
    from clavagen.examples import HBB_CDS

    def fake_auto(gene, tx, log=None):
        return {"gene": "HBB", "known_variants": [
            {"hgvs_c": "", "hgvs_p": "p.(Lys67Glu)", "classification": "Pathogenic", "source": "ClinVar"}]}

    monkeypatch.setattr(online, "auto_annotation", fake_auto)
    monkeypatch.setattr(online, "http", lambda *a, **k: (_ for _ in ()).throw(OSError("offline test")))
    ref = tmp_path / "ref.fa"
    ref.write_text(">HBB_test (HBB)\n" + HBB_CDS + "\n")
    query = HBB_CDS[:198] + "G" + HBB_CDS[199:]                       # c.199A>G p.(Lys67Glu)

    def met(exclude):
        p = pl_mod.Pipeline([f"q={query}"], [str(ref)], gene="HBB", online=True, use_clinvar=False,
                            auto_cache={}, log=lambda *_: None, exclude_known=exclude)
        v = p.run()["samples"][0]["variants"][0]
        return {c["code"] for c in v["criteria"] if c["status"] == "met"}

    assert "PS1" in met({"c.199A>G"})                       # DNA-only exclusion misses the protein-only record
    assert "PS1" not in met({"c.199A>G", "p.Lys67Glu"})     # as the benchmark now does: record left out


def test_old_transcript_version_still_gets_population_data(monkeypatch, tmp_path):
    """ERepo rows on an older RefSeq version (e.g. NM_000257.3) used to get no gnomAD data, because
    Ensembl VEP only resolves current versions: common benign variants then missed BA1."""
    from test_online import _genbank
    import urllib.parse

    def fake(url, data=None, **k):
        if "efetch.fcgi" in url:
            acc = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["id"][0]
            return _genbank().replace("NM_000518.5", "NM_000518.4") if acc == "NM_000518.4" else _genbank()
        if "rest.ensembl.org/vep" in url and any("NM_000518.4" in h for h in data["hgvs_notations"]):
            raise online.OnlineError("rest.ensembl.org: HTTP 400 (unknown transcript version)")
        return fake_http(url, data, **k)
    monkeypatch.setattr(online, "http", fake)
    online._dead_hosts.clear()
    old = EREPO.replace("NM_000518.5", "NM_000518.4")
    recs, _ = load_erepo(old)
    s = run_benchmark([r for r in recs if r["hgvs_c"] == "c.9T>C"], str(tmp_path), log=lambda *_: None)
    assert s["analysed"] == 1 and s["exact_agreement"] == 1.0          # BA1 -> Benign, as the experts said
    tsv = open(os.path.join(tmp_path, "benchmark_results.tsv")).read()
    assert "not the current RefSeq version" in tsv and "0.25" in tsv


# ---------------------------------------------------------------- publication outputs
def test_publication_outputs_and_resume(mocked, tmp_path):
    """Full run writes per-variant JSONL, statistics with CIs, manuscript tables and methods text;
    an interrupted run resumes without re-analysing finished variants."""
    import csv
    recs, stats = load_erepo(EREPO)
    meta = {"erepo_filter": stats, "erepo_file": "erepo.tsv", "erepo_sha256": "ab" * 32}
    s = run_benchmark(recs, str(tmp_path), log=lambda *_: None, meta=meta, figures=False)
    st = s["statistics"]
    assert st["three_tier"]["n"] == 3 and len(st["three_tier"]["ci95"]) == 2
    assert set(st["by_variant_type"]) >= {"missense", "nonsense"} and "Hemoglobinopathy VCEP" in st["by_panel"]
    assert s["run"]["erepo_filter"]["used"] == 3 and s["run"]["selected"] == 3 and s["run"]["finished"]
    for f in ("tables/table1_overall.csv", "tables/table2_by_panel.csv", "tables/table3_criteria.csv",
              "tables/table4_disagreement.csv", "tables/confusion_matrix.csv", "methods_and_results.md"):
        assert (tmp_path / f).exists(), f
    t1 = list(csv.reader(open(tmp_path / "tables/table1_overall.csv")))
    assert t1[1][0].startswith("3-tier") and "(" in t1[1][2]
    assert "SHA-256 abababab" in open(tmp_path / "methods_and_results.md").read()
    lines = open(tmp_path / "benchmark_results.jsonl").read().splitlines()
    assert len(lines) == 3
    # simulate an interruption after the first variant, then resume
    open(tmp_path / "benchmark_results.jsonl", "w").write(lines[0] + "\n" + lines[1][:20])
    seen = []
    from clavagen import pipeline
    orig = pipeline.Pipeline.run
    monkey = pytest.MonkeyPatch()
    monkey.setattr(pipeline.Pipeline, "run", lambda self: seen.append(1) or orig(self))
    s2 = run_benchmark(recs, str(tmp_path), log=lambda *_: None, resume=True, figures=False)
    monkey.undo()
    assert len(seen) == 2 and s2["analysed"] == 3 and s2["run"]["resumed_from"] == 1


def test_wilson_and_kappa():
    from clavagen.benchmark_stats import kappa, wilson
    p, lo, hi = wilson(81, 263)                     # reference values: 0.3080 (0.2553-0.3662)
    assert (round(lo, 3), round(hi, 3)) == (0.255, 0.366)
    assert wilson(0, 10)[1] == 0.0 and wilson(10, 10)[2] == 1.0
    pairs = [("a", "a")] * 20 + [("a", "b")] * 5 + [("b", "a")] * 10 + [("b", "b")] * 15
    assert round(kappa(pairs, ["a", "b"]), 3) == 0.4  # textbook example: po 0.70, pe 0.50
    perfect = [("P", "P"), ("V", "V"), ("B", "B")] * 4
    assert kappa(perfect, ["P", "V", "B"]) == 1.0
    assert kappa(perfect, ["P", "V", "B"], weighted=True) == 1.0


def test_figures_written(mocked, tmp_path):
    pytest.importorskip("matplotlib")
    recs, stats = load_erepo(EREPO)
    run_benchmark(recs, str(tmp_path), log=lambda *_: None, meta={"erepo_filter": stats})
    figs = sorted(p.name for p in (tmp_path / "figures").iterdir())
    for name in ("Fig1_variant_selection", "Fig2_confusion_matrix", "Fig3_agreement", "Fig5_criteria"):
        assert {f"{name}.png", f"{name}.pdf", f"{name}.svg"} <= set(figs), name
    assert "figure_legends.md" in figs
    from PIL import Image
    im = Image.open(tmp_path / "figures" / "Fig2_confusion_matrix.png")
    assert round(im.info["dpi"][0]) == 600


def test_benchmark_without_matplotlib(mocked, tmp_path, monkeypatch):
    """matplotlib is optional: results, tables and text are written, figures are skipped with a hint."""
    import sys
    monkeypatch.setitem(sys.modules, "matplotlib", None)
    msgs = []
    recs, _ = load_erepo(EREPO)
    s = run_benchmark(recs, str(tmp_path), log=msgs.append)
    assert s["analysed"] == 3 and (tmp_path / "tables" / "table1_overall.csv").exists()
    assert not (tmp_path / "figures").exists() or not any((tmp_path / "figures").iterdir())
    assert any("pip install matplotlib" in m for m in msgs)


def test_vep_falls_back_to_original_version(monkeypatch, tmp_path):
    """NCBI can be ahead of VEP (NM_000249.5 at NCBI, only .4 in VEP): retry the original version."""
    from test_online import _genbank
    import urllib.parse

    def fake(url, data=None, **k):
        if "efetch.fcgi" in url:
            acc = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["id"][0]
            return _genbank().replace("NM_000518.5", "NM_000518.6") if acc == "NM_000518" else _genbank()
        if "rest.ensembl.org/vep" in url and any("NM_000518.6" in h for h in data.get("hgvs_notations", [])):
            raise online.OnlineError("rest.ensembl.org: HTTP 400 (Could not get a Transcript object)")
        return fake_http(url, data, **k)
    monkeypatch.setattr(online, "http", fake)
    online._dead_hosts.clear()
    recs, _ = load_erepo(EREPO)
    run_benchmark([r for r in recs if r["hgvs_c"] == "c.9T>C"], str(tmp_path), log=lambda *_: None, figures=False)
    import csv
    row = next(csv.DictReader(open(tmp_path / "benchmark_results.tsv"), delimiter="\t"))
    assert row["vep_resolved"] == "True" and row["lookup_transcript"] == "NM_000518.5"
    assert row["tool_class"] == "Benign"                    # BA1 from the AF found on .5


def test_network_failure_is_retried_on_resume(monkeypatch, tmp_path):
    state = {"down": True}

    def fake(url, data=None, **k):
        if state["down"] and "rest.ensembl.org/vep" in url:
            raise online.OnlineError("rest.ensembl.org: no internet connection (server name could not be resolved)")
        return fake_http(url, data, **k)
    monkeypatch.setattr(online, "http", fake)
    online._dead_hosts.clear()
    recs, _ = load_erepo(EREPO)
    s = run_benchmark(recs, str(tmp_path), log=lambda *_: None, figures=False)
    assert s["errors"] == 3 and s["analysed"] == 0          # not silently classified without gnomAD
    state["down"] = False
    online._dead_hosts.clear()
    s = run_benchmark(recs, str(tmp_path), log=lambda *_: None, figures=False, resume=True)
    assert s["errors"] == 0 and s["analysed"] == 3


def test_shifted_indel_is_recovered(mocked, tmp_path):
    """ClinGen sometimes gives an indel at its 5' position (c.904delT); CLAVAGEN reports the 3'-shifted HGVS."""
    from clavagen.benchmark import apply_c
    from test_online import _genbank
    from clavagen.online import parse_genbank
    gb = parse_genbank(_genbank())
    cds = gb["seq"][gb["cds"][0] - 1:gb["cds"][1]]
    i = next(i for i in range(30, len(cds) - 2) if cds[i] == cds[i + 1] != cds[i - 1])   # a 2-base run
    five, three = f"c.{i + 1}del", f"c.{i + 2}del"
    assert apply_c(gb["seq"], gb["cds"], five)[0] == apply_c(gb["seq"], gb["cds"], three)[0]
    row = "\t".join(["x", "1", "CA1", f"NM_000518.5:{five}", "HBB", "beta thalassemia", "", "AR", "Pathogenic",
                     "PVS1, PM2_Supporting", "", "Hemoglobinopathy VCEP", "false"])
    recs, _ = load_erepo(EREPO.split("\n")[0] + "\n" + row + "\n")
    s = run_benchmark(recs, str(tmp_path), log=lambda *_: None, figures=False)
    assert s["analysed"] == 1 and s["not_recovered"] == 0


def test_not_testable_is_not_an_error():
    from clavagen.benchmark import summarise
    rows = [{"status": "error", "note": "c.-1170C>T outside transcript"},
            {"status": "error", "note": "reference base mismatch for c.*2050T>A (transcript has A)"},
            {"status": "error", "note": "network: no internet connection"}]
    s = summarise(rows)
    assert s["not_testable"] == 2 and s["errors"] == 1
