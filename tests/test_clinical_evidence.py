"""Guided clinical evidence and gene rule profiles (CLAVAGEN 1.0.2)."""

import json
import os
import subprocess

import pytest

from clavagen.acmg import Criterion, apply_gene_rules, combine
from clavagen.evidence import build
from clavagen.examples import write_examples
from clavagen.pipeline import Pipeline


def _one(d, rules=None):
    out = build(d, rules)
    assert len(out) == 1, out
    return out[0][:2]


# ---------------------------------------------------------------- conversion rules
def test_segregation_meioses():
    assert _one({"segregation": {"informative_meioses": 2}}) == ("PP1", None)
    assert _one({"segregation": {"informative_meioses": 3}}) == ("PP1", "Supporting")
    assert _one({"segregation": {"informative_meioses": 6}}) == ("PP1", "Moderate")
    assert _one({"segregation": {"informative_meioses": 7}}) == ("PP1", "Strong")
    assert _one({"segregation": {"lod": 0.9}}) == ("PP1", "Supporting")          # 3 x log10(2)
    # a panel's own thresholds (e.g. Jarvik & Browning single-family 1/32 = 5 meioses for Strong)
    rules = {"source": "panel X", "pp1_meioses": {"Supporting": 3, "Moderate": 4, "Strong": 5}}
    assert _one({"segregation": {"informative_meioses": 5}}, rules) == ("PP1", "Strong")


def test_pm3_points():
    p = {"partner": "Pathogenic", "phase": "confirmed"}
    assert _one({"in_trans": [p]}) == ("PM3", "Moderate")                         # 1.0
    assert _one({"in_trans": [p, p]}) == ("PM3", "Strong")                        # 2.0
    assert _one({"in_trans": [p] * 4}) == ("PM3", "VeryStrong")                   # 4.0
    assert _one({"in_trans": [{"partner": "Pathogenic", "phase": "unknown"}]}) == ("PM3", "Supporting")   # 0.5
    assert _one({"in_trans": [{"partner": "Likely pathogenic", "phase": "unknown"}]}) == ("PM3", None)     # 0.25
    assert _one({"homozygous_probands": 5}) == ("PM3", "Moderate")                # capped at 1.0


def test_functional_oddspath():
    assert _one({"functional": {"result": "damaging", "oddspath": 20}}) == ("PS3", "Strong")
    assert _one({"functional": {"result": "damaging", "oddspath": 5}}) == ("PS3", "Moderate")
    assert _one({"functional": {"result": "damaging", "oddspath": 2.5}}) == ("PS3", "Supporting")
    assert _one({"functional": {"result": "damaging", "oddspath": 1.5}}) == ("PS3", None)
    assert _one({"functional": {"result": "normal", "oddspath": 0.04}}) == ("BS3", "Strong")
    assert _one({"functional": {"result": "normal", "oddspath": 0.3}}) == ("BS3", "Supporting")
    assert _one({"functional": {"result": "damaging", "strength": "Moderate"}}) == ("PS3", "Moderate")


def test_de_novo_points():
    assert _one({"de_novo": [{"confirmed": True, "phenotype": "specific"}]}) == ("PS2", "Strong")      # 2
    assert _one({"de_novo": [{"confirmed": False, "phenotype": "specific"}]}) == ("PM6", "Moderate")   # 1
    assert _one({"de_novo": [{"confirmed": True, "phenotype": "specific"}] * 2}) == ("PS2", "VeryStrong")
    # low-specificity phenotypes count for at most 1 point in total
    assert _one({"de_novo": [{"confirmed": True, "phenotype": "heterogeneous"}] * 6}) == ("PS2", "Moderate")


def test_other_codes():
    codes = dict((c, s) for c, s, _ in build({"phenotype_specific": True, "non_segregation": True,
                                               "benign_phase": True, "alternative_cause": True,
                                               "healthy_adults": {"count": 2},
                                               "case_control": {"odds_ratio": 6, "ci_lower": 1.5}}))
    assert codes == {"PP4": "Supporting", "BS4": "Strong", "BP2": "Supporting", "BP5": "Supporting",
                     "BS2": "Strong", "PS4": "Strong"}
    assert _one({"case_control": {"odds_ratio": 3, "ci_lower": 1.2}}) == ("PS4", None)


# ---------------------------------------------------------------- gene rules
def test_gene_rules_disable_fix_cap():
    C = [Criterion("PP2", "met", "Supporting"), Criterion("PM2", "met", "Supporting"),
         Criterion("PP3", "met", "Strong"), Criterion("PVS1", "met", "VeryStrong")]
    rules = {"source": "Test VCEP v1", "disabled": ["PP2"], "strength": {"PM2": "Moderate"},
             "max_strength": {"PP3": "Supporting"}}
    out = {c.code: (c.status, c.strength) for c in apply_gene_rules(C, {"rules": rules})}
    assert out == {"PP2": ("not_met", ""), "PM2": ("met", "Moderate"), "PP3": ("met", "Supporting"),
                   "PVS1": ("met", "VeryStrong")}


# ---------------------------------------------------------------- end to end
@pytest.fixture
def demo(tmp_path):
    return write_examples(str(tmp_path))


def _criteria(res, hgvs):
    v = next(v for s in res["samples"] for v in s["variants"] if v["hgvs_c"] == hgvs)
    return v, {c["code"]: c for c in v["criteria"] if c["status"] == "met"}


def test_pipeline_uses_clinical_evidence(demo):
    base = Pipeline([demo["dna_queries"]], [demo["reference_dna"]], annotation=demo["annotation"]).run()
    v0, c0 = _criteria(base, "c.20A>T")
    ce = {"c.20A>T": {"functional": {"result": "damaging", "oddspath": 25},
                      "in_trans": [{"partner": "Pathogenic", "phase": "confirmed"}] * 2}}
    res = Pipeline([demo["dna_queries"]], [demo["reference_dna"]], annotation=demo["annotation"],
                   clinical_evidence=ce).run()
    v, c = _criteria(res, "c.20A>T")
    assert c["PS3"]["strength"] == "Strong" and c["PS3"]["source"] == "clinical"
    assert c["PM3"]["strength"] == "Strong" and "2 PM3 points" in c["PM3"]["rationale"]
    assert v["points"] > v0["points"]
    # other variants are not affected by a variant-specific entry
    _, c118 = _criteria(res, "c.118C>T")
    assert "PS3" not in c118


def test_pipeline_uses_gene_rules(demo):
    ann = json.load(open(demo["annotation"]))
    ann["rules"] = {"source": "Test panel", "disabled": ["PM5"], "max_strength": {"PP3": "Supporting"}}
    res = Pipeline([demo["dna_queries"]], [demo["reference_dna"]], annotation=ann).run()
    v, c = _criteria(res, "c.20A>T")
    assert "PM5" not in c and c["PP3"]["strength"] == "Supporting"
    pm5 = next(x for x in v["criteria"] if x["code"] == "PM5")
    assert "Not used for this gene (Test panel)" in pm5["rationale"]


def test_cli_rules_and_evidence_files(demo, tmp_path):
    from clavagen.cli import _apply_rules
    rp = tmp_path / "rules.json"
    rp.write_text(json.dumps({"HBB": {"source": "X", "thresholds": {"pm2": 0.0}, "disabled": ["PP3"],
                                       "_help": "ignored"}}))
    ann = _apply_rules({"thresholds": {"ba1": 0.05}}, str(rp))
    assert ann["references"]["HBB"]["rules"] == {"source": "X", "disabled": ["PP3"]}
    assert ann["references"]["HBB"]["thresholds"] == {"pm2": 0.0}
    single = tmp_path / "one.json"
    single.write_text(json.dumps({"source": "Y", "max_strength": {"PP3": "Supporting"}, "thresholds": {"bs1": 0.001}}))
    ann = _apply_rules({"thresholds": {"ba1": 0.05}}, str(single))
    assert ann["rules"]["max_strength"] == {"PP3": "Supporting"} and ann["thresholds"] == {"ba1": 0.05, "bs1": 0.001}


def test_cli_run_with_clinical_evidence(demo, tmp_path):
    from clavagen.cli import main
    ce = tmp_path / "ce.json"
    ce.write_text(json.dumps({"c.20A>T": {"phenotype_specific": True}}))
    main(["run", "-q", demo["dna_queries"], "-r", demo["reference_dna"], "-a", demo["annotation"], "--offline",
          "--clinical-evidence", str(ce), "-o", str(tmp_path / "out"), "--formats", "json"])
    js = next(f for f in os.listdir(tmp_path / "out") if f.endswith(".json"))
    _, c = _criteria(json.load(open(tmp_path / "out" / js)), "c.20A>T")
    assert c["PP4"]["source"] == "clinical"


def test_patient_table_clinical_column(tmp_path, demo):
    from clavagen.cohort import read_patient_table
    q = open(demo["dna_queries"]).read().split(">")[1]
    rows = [{"id": "P1", "query": ">" + q, "clinical_evidence": json.dumps({"*": {"phenotype_specific": True}})}]
    pts = read_patient_table(json.dumps(rows))
    assert pts[0]["clinical_evidence"] == {"*": {"phenotype_specific": True}}
    with pytest.raises(ValueError):
        read_patient_table(json.dumps([{"id": "P2", "query": ">" + q, "clinical_evidence": "{bad"}]))


def test_app_form_builds_evidence_json():
    """The app window's guided form (JavaScript) turns the fields into the JSON the engine reads, and back."""
    import re
    import shutil
    from clavagen.webui import PAGE
    if not shutil.which("node"):
        pytest.skip("node not installed")
    js = "\n".join(re.findall(r"<script>(.*?)</script>", PAGE, re.S))
    funcs = "\n".join(re.findall(r"(const CEN=.*?\nfunction ceSummary.*?\n)(?=\$|function |const |let )", js, re.S))
    assert "function readCE" in funcs
    harness = """
const el={};const $=id=>el[id]||(el[id]={value:"",checked:false,open:false});
""" + funcs + """
$("ce-key").value="c.20A>T";$("ce-meioses").value="5";$("ce-tp-c").value="2";$("ce-func").value="damaging";
$("ce-oddspath").value="20";$("ce-dn-a").value="1";$("ce-dn-ph").value="specific";$("ce-pp4").checked=true;$("ce-apstr").value="Supporting";
const ce=readCE();fillCE(ce);const again=readCE();
console.log(JSON.stringify({ce,again,sum:ceSummary(ce)}));
"""
    out = subprocess.run(["node", "-e", harness], capture_output=True, text=True, check=True).stdout
    d = json.loads(out)
    e = d["ce"]["c.20A>T"]
    assert e["segregation"]["informative_meioses"] == 5 and len(e["in_trans"]) == 2
    assert e["functional"] == {"result": "damaging", "oddspath": 20} and e["phenotype_specific"] is True
    assert e["de_novo"] == [{"confirmed": False, "phenotype": "specific"}]
    assert d["again"] == d["ce"]                                    # edit -> save keeps everything
    codes = {c for c, s, _ in build(e)}
    assert codes == {"PP1", "PM3", "PS3", "PM6", "PP4"}


def test_app_payload_with_rules_and_clinical_evidence(tmp_path):
    """What the app window sends: patients with guided evidence + a gene rule profile box."""
    from clavagen.webui import _example_payload, run_payload
    d = _example_payload()
    d["online"] = "0"
    d["rules"] = json.dumps({"source": "Test panel", "max_strength": {"PP3": "Supporting"}})
    pt = next(p for p in d["patients"] if not p.get("evidence"))   # reviewer codes would take priority
    pt["clinical_evidence"] = {"*": {"phenotype_specific": True,
                                     "functional": {"result": "damaging", "oddspath": 30}}}
    summ = run_payload(d, str(tmp_path))
    found, pp3 = [], []
    for root, _, files in os.walk(tmp_path):
        for f in files:
            if f == "clinical_summary.json" and os.path.basename(root) == pt["id"]:
                try:
                    js = json.load(open(os.path.join(root, f)))
                except Exception:
                    continue
                for s in js.get("samples", []) if isinstance(js, dict) else []:
                    for v in s.get("variants", []):
                        found += [(c["code"], c["source"]) for c in v["criteria"] if c["status"] == "met"]
                        pp3 += [c["strength"] for c in v["criteria"] if c["status"] == "met" and c["code"] == "PP3"]
    assert ("PS3", "clinical") in found and ("PP4", "clinical") in found
    assert set(pp3) <= {"Supporting"}                              # rule profile cap applied
    assert summ
