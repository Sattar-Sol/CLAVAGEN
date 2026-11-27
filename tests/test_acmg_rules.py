"""ACMG rule fixes found in the full ClinGen benchmark (CLAVAGEN 1.0.2)."""

from clavagen.acmg import evaluate
from clavagen.consequence import _aa_repeat


def _v(cons, hgvs="c.100A>T", **kw):
    v = {"consequence": cons, "hgvs_c": hgvs, "hgvs_p": "p.(Lys34Ter)", "mol_type": "dna", "aa_pos": 34,
         "ptc_pos": 34, "domains": [], "splice_region": False, "kind": "SNV"}
    v.update(kw)
    return v


def _ctx(revel=None, spliceai=None, af=None, lof=True, hgvs="c.100A>T", **kw):
    ins = {k: x for k, x in (("REVEL", revel), ("SpliceAI", spliceai)) if x is not None}
    ctx = {"lof_mechanism": lof, "insilico": {hgvs: ins} if ins else {}, "population_dataset": "gnomAD",
           "population": {hgvs: {"af": af}} if af is not None else {}, "cds_length": 3000, "moi": "AD"}
    ctx.update(kw)
    return ctx


def _met(C):
    return {c.code: c.strength for c in C if c.status == "met"}


def test_no_splice_bp4_on_null_variants():
    """SpliceAI 0 says nothing about a nonsense or frameshift variant (BP4 was applied to 1,136 of them)."""
    for cons in ("stop_gained", "frameshift_variant", "start_lost", "inframe_deletion", "5_prime_UTR_variant"):
        met = _met(evaluate(_v(cons), _ctx(spliceai=0.0)))
        assert "BP4" not in met, cons
    assert "BP4" in _met(evaluate(_v("synonymous_variant"), _ctx(spliceai=0.0)))


def test_pp3_not_combined_with_pvs1():
    met = _met(evaluate(_v("stop_gained"), _ctx(spliceai=0.6)))
    assert "PVS1" in met and "PP3" not in met


def test_missense_bp4_needs_no_splice_impact():
    """REVEL benign but SpliceAI high: splicing PP3, not a contradictory BP4 + PP3 pair."""
    met = _met(evaluate(_v("missense_variant"), _ctx(revel=0.12, spliceai=0.85)))
    assert "BP4" not in met and met.get("PP3") == "Supporting"
    met = _met(evaluate(_v("missense_variant"), _ctx(revel=0.12, spliceai=0.15)))
    assert "BP4" not in met and "PP3" not in met            # 0.1-0.2: neither
    met = _met(evaluate(_v("missense_variant"), _ctx(revel=0.12, spliceai=0.02)))
    assert met.get("BP4") == "Moderate"


def test_bp3_needs_a_protein_repeat():
    # p.Phe299del in Phe-Phe: a nucleotide poly-T, but no protein repeat -> PM4, not BP3
    v = _v("inframe_deletion", "c.895_897del", repeat_unit="T", repeat_copies=5, kind="deletion",
           **_aa_repeat("MKFFLA", 2, "F"))
    met = _met(evaluate(v, _ctx(hgvs="c.895_897del")))
    assert "PM4" in met and "BP3" not in met
    # p.His20del in a run of 8 His -> BP3
    v = _v("inframe_deletion", "c.60_62del", kind="deletion", **_aa_repeat("MAA" + "H" * 8 + "GG", 5, "H"))
    met = _met(evaluate(v, _ctx(hgvs="c.60_62del")))
    assert "BP3" in met and "PM4" not in met


def test_aa_repeat_units():
    assert _aa_repeat("MSAGAGAGAGT", 3, "AG")["aa_repeat_copies"] == 4
    assert _aa_repeat("MKKLL", 1, "K")["aa_repeat_len"] == 2


def test_unresolved_variant_gets_no_pm2():
    """Another variant's successful lookup must not make this one look 'absent from gnomAD'."""
    ctx = _ctx(af=0.3, hgvs="c.5G>A")
    ctx["population_unresolved"] = {"c.100A>T"}
    C = evaluate(_v("missense_variant"), ctx)
    assert "PM2" not in _met(C)
    assert any(c.code == "PM2" and c.status == "not_evaluated" for c in C)


def test_lof_mechanism_counts_all_transcripts():
    """ClinGen variants on a non-MANE transcript (e.g. HNF1A NM_001306179.2) got no PVS1: only ClinVar records
    on that exact transcript were counted, and those are rarely truncating."""
    from clavagen.online import derive_gene_evidence, recalc_gene_evidence
    ann = {"known_variants": [{"hgvs_c": "c.10A>G", "hgvs_p": "p.(Lys4Arg)"}], "_gene_truncating_all": 5}
    derive_gene_evidence(ann)
    assert ann["lof_mechanism"] is True
    # leave-one-out: the tested variant may be one of them under another transcript's c. -> count minus one
    ctx = {"known_variants": [], "_gene_truncating_all": 3, "lof_mechanism": True, "_lof_derived": "clinvar"}
    recalc_gene_evidence(ctx)
    assert not ctx["lof_mechanism"]
    ctx = {"known_variants": [], "_gene_truncating_all": 4, "lof_mechanism": True, "_lof_derived": "clinvar"}
    recalc_gene_evidence(ctx)
    assert ctx["lof_mechanism"] is True


def test_auto_annotation_keeps_gene_level_truncating_count(monkeypatch):
    from clavagen import online
    known = [{"hgvs_c": f"c.{i}C>T", "hgvs_p": f"p.(Gln{i}Ter)", "transcript": "NM_000545.8"} for i in range(1, 6)]
    monkeypatch.setattr(online, "fetch_clinvar", lambda g, log=None: known)
    ann = online.auto_annotation("HNF1A", transcript="NM_001306179.2", sources=("clinvar",))
    assert ann["known_variants"] == [] and ann["_gene_truncating_all"] == 5 and ann["lof_mechanism"] is True
