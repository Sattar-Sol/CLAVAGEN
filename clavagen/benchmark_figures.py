"""Publication figures for the ClinGen benchmark (needs matplotlib).

Each figure is written three ways in ``figures/``:
  * ``.png`` at 600 dpi (journal raster requirement for line art is usually 600-1000 dpi),
  * ``.pdf`` and ``.svg`` (vector: scale to any size, edit text in Illustrator/Inkscape).

Layout follows common journal specifications: single column 89 mm, double column 183 mm,
sans-serif text 7-8 pt, colour-blind-safe palette (Okabe & Ito), no titles inside the figure
(the caption carries them; draft captions are written to ``figure_legends.md``).

Regenerate at any time from a finished run:  clavagen benchmark-figures <benchmark folder>
"""

from __future__ import annotations

import os

from .benchmark_stats import CLASSES, TIERS, prop, tier

MM = 1 / 25.4
SINGLE, DOUBLE = 89 * MM, 183 * MM
# Okabe & Ito colour-blind-safe palette
BLUE, ORANGE, GREEN, VERMIL, SKY, PURPLE, YELLOW, GREY = ("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9",
                                                           "#CC79A7", "#F0E442", "#8C8C8C")
TIER_COL = {"P/LP": VERMIL, "VUS": GREY, "B/LB": BLUE}
CAUSE_COL = {"missing evidence": GREY, "criterion assignment": ORANGE, "evidence strength": PURPLE}
SHORT = {"Pathogenic": "P", "Likely pathogenic": "LP", "Uncertain significance": "VUS", "Likely benign": "LB",
         "Benign": "B"}


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",   # text stays editable in vector files
        "savefig.dpi": 600, "figure.dpi": 150, "axes.titleweight": "bold", "axes.titlelocation": "left"})
    return plt


def _save(fig, out, name, made):
    for ext in ("png", "pdf", "svg"):
        p = os.path.join(out, f"{name}.{ext}")
        fig.savefig(p, bbox_inches="tight", pad_inches=0.03, dpi=600 if ext == "png" else None)
        made.append(p)
    import matplotlib.pyplot as plt
    plt.close(fig)


def _panel(ax, letter):
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")


def _pct(d):
    return 100 * d["p"] if d and d.get("p") is not None else 0.0


def _err(d):
    if not d or d.get("p") is None:
        return [[0], [0]]
    return [[100 * (d["p"] - d["ci95"][0])], [100 * (d["ci95"][1] - d["p"])]]


# ------------------------------------------------------------------ figure 1: variant selection
def fig_selection(plt, summary, out, made):
    run = summary.get("run") or {}
    f = run.get("erepo_filter") or {}
    if not f:
        return
    steps = [("Evidence Repository rows", f.get("total", 0))]
    left = f.get("total", 0)
    for lab, k in (("retracted", "retracted"), ("no 5-tier classification", "no_class"),
                   ("no NM_ c. description", "no_transcript_hgvs"), ("intronic / complex", "intronic_or_complex")):
        left -= f.get(k, 0)
        steps.append((f"minus {lab} ({f.get(k, 0):,})", left))
    if run.get("genes") not in (None, "all") or run.get("max"):
        steps.append(("selected for this run", run.get("selected", 0)))
    sel = steps[-1][1]
    for lab, k in (("not on transcript record", "not_testable"), ("not recovered as one call", "not_recovered"),
                   ("failed", "errors")):
        if summary.get(k):
            sel -= summary[k]
            steps.append((f"minus {lab} ({summary[k]:,})", sel))
    steps.append(("analysed", summary["analysed"]))
    fig, ax = plt.subplots(figsize=(SINGLE, 0.28 * len(steps) + 0.35))
    ys = range(len(steps))[::-1]
    mx = max(v for _, v in steps) or 1
    ax.barh(list(ys), [v for _, v in steps], color=[GREY] * (len(steps) - 1) + [BLUE], height=0.62)
    for y, (lab, v) in zip(ys, steps):
        ax.text(v + mx * 0.01, y, f"{v:,}", va="center", fontsize=6.5)
    ax.set_yticks(list(ys), [lab for lab, _ in steps])
    ax.set_xlim(0, mx * 1.15)
    ax.set_xlabel("Variants")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    _save(fig, out, "Fig1_variant_selection", made)


# ------------------------------------------------------------------ figure 2: confusion matrix
def fig_confusion(plt, summary, out, made):
    import numpy as np
    conf = summary["confusion"]
    m = np.array([[conf[e][t] for t in CLASSES] for e in CLASSES], float)
    rows = m.sum(axis=1, keepdims=True)
    frac = np.divide(m, rows, out=np.zeros_like(m), where=rows > 0)
    fig, ax = plt.subplots(figsize=(SINGLE, SINGLE * 0.86))
    im = ax.imshow(frac, cmap="Blues", vmin=0, vmax=1)
    for i in range(5):
        for j in range(5):
            c = int(m[i, j])
            ax.text(j, i, f"{c}\n{100 * frac[i, j]:.0f}%" if rows[i, 0] else "0", ha="center", va="center",
                    fontsize=6, color="white" if frac[i, j] > 0.55 else "black",
                    fontweight="bold" if i == j else "normal")
    labels = [SHORT[c] for c in CLASSES]
    ax.set_xticks(range(5), labels)
    ax.set_yticks(range(5), [f"{lab} (n={int(rows[i, 0])})" for i, lab in enumerate(labels)])
    ax.set_xlabel("CLAVAGEN classification")
    ax.set_ylabel("Expert panel classification")
    for x in (1.5, 2.5):   # tier boundaries
        ax.axhline(x, color="black", lw=0.5)
        ax.axvline(x, color="black", lw=0.5)
    ax.spines[:].set_visible(False)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Share of expert row", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_linewidth(0.5)
    _save(fig, out, "Fig2_confusion_matrix", made)


# ------------------------------------------------------------------ figure 3: agreement with CIs
def fig_agreement(plt, summary, out, made):
    st = summary["statistics"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 2.55), gridspec_kw={"width_ratios": [1.15, 1]})
    items = [("3-tier agreement", st["three_tier"], BLUE),
             ("3-tier, expert class from\ncomputable criteria", st["computable_three_tier"], SKY),
             ("Exact 5-class agreement", st["exact"], BLUE),
             ("Exact, computable criteria", st["computable_exact"], SKY),
             *([("3-tier, experts' clinical\nevidence entered (simulation)", st["with_clinical_evidence"]["three_tier"],
                 GREEN)] if st.get("with_clinical_evidence") else []),
             ("Critical discordance\n(P/LP vs B/LB)", st["critical_discordance"], VERMIL)]
    ys = list(range(len(items)))[::-1]
    for y, (lab, d, col) in zip(ys, items):
        a.barh(y, _pct(d), color=col, height=0.6)
        a.errorbar(_pct(d), y, xerr=_err(d), fmt="none", ecolor="black", elinewidth=0.6, capsize=1.8)
        a.text(min(_pct(d) + (d["ci95"][1] - d["p"]) * 100 + 2, 101), y, f"{_pct(d):.1f}%", va="center", fontsize=6.5)
    a.set_yticks(ys, [lab for lab, _, _ in items])
    a.set_xlim(0, 112)
    a.set_xticks([0, 25, 50, 75, 100])
    a.set_xlabel("% of variants (95% CI)")
    a.tick_params(axis="y", length=0)
    k3, k5 = st["kappa_three_tier"], st["weighted_kappa_five_class"]
    if k3.get("kappa") is not None and k3["ci95"][0] is not None and k5.get("kappa") is not None \
            and k5["ci95"][0] is not None:
        a.text(1.0, -0.3, f"Cohen's κ (3-tier) {k3['kappa']:.2f} [{k3['ci95'][0]:.2f}, {k3['ci95'][1]:.2f}];  "
               f"weighted κ (5-class) {k5['kappa']:.2f} [{k5['ci95'][0]:.2f}, {k5['ci95'][1]:.2f}]",
               transform=a.transAxes, ha="right", fontsize=6, color="#444")
    _panel(a, "a")
    bt = st["by_expert_tier"]
    xs = [t for t in TIERS if t in bt]
    vals = [bt[t]["three_tier"] for t in xs]
    b.bar(range(len(xs)), [_pct(d) for d in vals], color=[TIER_COL[t] for t in xs], width=0.6)
    b.errorbar(range(len(xs)), [_pct(d) for d in vals], yerr=[[e[0][0] for e in map(_err, vals)],
                                                              [e[1][0] for e in map(_err, vals)]],
               fmt="none", ecolor="black", elinewidth=0.6, capsize=2)
    for i, d in enumerate(vals):
        b.text(i, _pct(d) + (d["ci95"][1] - d["p"]) * 100 + 2.5, f"{_pct(d):.0f}%", ha="center", fontsize=6.5)
    b.set_xticks(range(len(xs)), [f"Expert {t}\n(n={bt[t]['n']})" for t in xs])
    b.set_ylim(0, 112)
    b.set_yticks([0, 25, 50, 75, 100])
    b.set_ylabel("CLAVAGEN in the same tier (%)")
    _panel(b, "b")
    fig.tight_layout(w_pad=2.5)
    _save(fig, out, "Fig3_agreement", made)


# ------------------------------------------------------------------ figure 4: forest plot by panel
def fig_panels(plt, summary, out, made, min_n=10):
    st = summary["statistics"]
    groups = [(g, d) for g, d in st["by_panel"].items() if d["n"] >= min_n and g != "unknown"]
    if len(groups) < 2:
        groups = [(g, d) for g, d in st["by_gene"].items() if d["n"] >= min_n] or list(st["by_gene"].items())
        what = "gene"
    else:
        what = "panel"
    if not groups:
        return
    groups.sort(key=lambda gd: gd[1]["three_tier"]["p"] or 0)
    small = sum(d["n"] for g, d in st["by_" + what].items() if d["n"] < min_n)
    h = max(1.6, 0.17 * len(groups) + 0.7)
    fig, ax = plt.subplots(figsize=(DOUBLE * 0.8, h))
    for y, (g, d) in enumerate(groups):
        p = _pct(d["three_tier"])
        ax.errorbar(p, y, xerr=_err(d["three_tier"]), fmt="s", ms=2.8 + min(3.5, d["n"] / 40), color=BLUE,
                    ecolor=BLUE, elinewidth=0.8, capsize=0)
        if d["critical"]["k"]:
            ax.text(103, y, f"{d['critical']['k']} critical", va="center", fontsize=5.5, color=VERMIL)
    overall = _pct(st["three_tier"])
    ax.axvline(overall, color=GREY, lw=0.7, ls="--")
    ax.text(overall, len(groups) - 0.3, f" overall {overall:.0f}%", fontsize=6, color="#555", va="bottom")
    ax.set_yticks(range(len(groups)), [f"{_short(g)} (n={d['n']})" for g, d in groups])
    ax.set_xlim(0, 115)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("3-tier agreement with expert classification, % (95% CI)")
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.7, len(groups) - 0.1)
    if small:
        ax.annotate(f"{what.capitalize()}s with fewer than {min_n} variants not shown ({small} variants).",
                    xy=(0, 0), xycoords="axes fraction", xytext=(0, -30), textcoords="offset points",
                    fontsize=5.5, color="#555", va="top")
    _save(fig, out, f"Fig4_agreement_by_{what}", made)


def _short(name):
    name = " ".join(name.replace("Variant Curation Expert Panel", "VCEP").replace("Expert Panel", "EP").split())
    return name if len(name) <= 42 else name[:40] + "…"


# ------------------------------------------------------------------ figure 5: criteria
def fig_criteria(plt, summary, out, made, min_obs=5):
    from .benchmark import COMPUTABLE
    crit = summary["criteria"]
    comp = [(c, d) for c, d in crit.items() if c in COMPUTABLE and d["tp"] + d["fp"] + d["fn"] >= min_obs]
    nonc = [(c, d["tp"] + d["fn"]) for c, d in crit.items() if c not in COMPUTABLE and d["tp"] + d["fn"] > 0]
    nonc.sort(key=lambda cd: -cd[1])
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, max(2.0, 0.2 * max(len(comp), len(nonc[:12])) + 0.8)),
                               gridspec_kw={"width_ratios": [1.3, 1]})
    comp.sort(key=lambda cd: cd[0])
    ys = list(range(len(comp)))[::-1]
    for y, (c, d) in zip(ys, comp):
        for val, n, col, off, mk in ((d["tp"], d["tp"] + d["fp"], ORANGE, 0.17, "o"), (d["tp"], d["tp"] + d["fn"], BLUE, -0.17, "D")):
            pr = prop(val, n)
            if pr["p"] is None:
                continue
            a.errorbar(_pct(pr), y + off, xerr=_err(pr), fmt=mk, ms=3, color=col, ecolor=col, elinewidth=0.7)
    a.set_yticks(ys, [f"{c} ({d['tp'] + d['fn']}/{d['tp'] + d['fp']})" for c, d in comp])
    a.set_xlim(-2, 102)
    a.set_xlabel("% (95% CI)")
    a.tick_params(axis="y", length=0)
    a.plot([], [], "o", color=ORANGE, ms=3, label="Precision (tool-applied also applied by experts)")
    a.plot([], [], "D", color=BLUE, ms=3, label="Recall (expert-applied also applied by tool)")
    a.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), frameon=False, ncol=1)
    a.set_title("Computable criteria (expert n / tool n)", fontsize=7)
    _panel(a, "a")
    nonc = nonc[:12]
    n_ok = summary["analysed"] or 1
    by = list(range(len(nonc)))[::-1]
    b.barh(by, [100 * n / n_ok for _, n in nonc], color=GREY, height=0.6)
    for y, (c, n) in zip(by, nonc):
        b.text(100 * n / n_ok + 0.8, y, f"{n}", va="center", fontsize=6)
    b.set_yticks(by, [c for c, _ in nonc])
    b.set_xlabel("% of variants where experts applied it")
    b.tick_params(axis="y", length=0)
    b.set_title("Criteria a sequence-based tool cannot apply", fontsize=7)
    _panel(b, "b")
    fig.tight_layout(w_pad=2.5)
    _save(fig, out, "Fig5_criteria", made)


# ------------------------------------------------------------------ figure 6: causes of disagreement
def fig_disagreement(plt, results, out, made):
    ok = [r for r in results if r["status"] == "ok" and r["tool_class"] != r["expert_class"]]
    if not ok:
        return
    causes = list(CAUSE_COL)
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 2.0), gridspec_kw={"width_ratios": [1, 1.1]})
    bottom = [0] * 3
    for cause in causes:
        vals = [sum(1 for r in ok if tier(r["expert_class"]) == t and r.get("disagreement") == cause) for t in TIERS]
        a.bar(range(3), vals, bottom=bottom, color=CAUSE_COL[cause], width=0.6, label=cause)
        bottom = [x + y for x, y in zip(bottom, vals)]
    for i, tot in enumerate(bottom):
        a.text(i, tot + max(bottom) * 0.02, str(tot), ha="center", fontsize=6.5)
    a.set_xticks(range(3), [f"Expert {t}" for t in TIERS])
    a.set_ylabel("Discordant variants (5-class)")
    from matplotlib.ticker import MaxNLocator
    a.yaxis.set_major_locator(MaxNLocator(integer=True))
    a.set_ylim(0, max(bottom) * 1.12 + 1)
    a.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, fontsize=6, handlelength=1,
             columnspacing=1.0)
    _panel(a, "a")
    from collections import Counter
    over = Counter((SHORT[r["expert_class"]], SHORT[r["tool_class"]]) for r in ok)
    top = over.most_common(10)[::-1]
    ys = range(len(top))
    i = {SHORT[c]: k for k, c in enumerate(CLASSES)}
    b.barh(list(ys), [n for _, n in top], color=[VERMIL if i[t] < i[e] else BLUE for (e, t), _ in top], height=0.6)
    for y, (_, n) in zip(ys, top):
        b.text(n + 0.3, y, str(n), va="center", fontsize=6)
    b.set_yticks(list(ys), [f"expert {e} → tool {t}" for (e, t), _ in top])
    b.tick_params(axis="y", length=0)
    b.set_xlabel("Variants")
    b.plot([], [], "s", color=VERMIL, label="tool more pathogenic")
    b.plot([], [], "s", color=BLUE, label="tool more benign")
    b.legend(frameon=False, loc="lower right", fontsize=6)
    _panel(b, "b")
    fig.tight_layout(w_pad=2.5)
    _save(fig, out, "Fig6_disagreement", made)


LEGENDS = """# Draft figure legends

Adjust wording and numbering to the journal. Numbers come from benchmark_summary.json.

**Figure 1. Selection of expert-classified variants.** Variants in the ClinGen Evidence Repository export ({total:,} rows) were excluded if retracted, without a 5-tier classification, without an NM_ transcript c. description, or intronic/complex; {analysed:,} variants were analysed.

**Figure 2. Classification by CLAVAGEN versus ClinGen expert panels.** Counts and row percentages (share of each expert class). Lines mark the 3-tier boundaries (P/LP, VUS, B/LB); bold cells on the diagonal are exact agreement. n = {analysed:,} variants.

**Figure 3. Agreement with expert classifications.** (a) Agreement and critical discordance with Wilson 95% confidence intervals; "computable criteria" = expert classification recalculated from the criteria a sequence-based tool can evaluate; "clinical evidence entered (simulation)" = CLAVAGEN's automated criteria plus the experts' clinical, functional, segregation and case criteria, as if entered with the guided evidence form. Cohen's kappa and linearly weighted kappa with bootstrap 95% intervals. (b) 3-tier agreement by expert tier.

**Figure 4. Agreement by expert panel (or gene).** 3-tier agreement with Wilson 95% confidence intervals; marker size reflects the number of variants; dashed line, overall agreement. Groups with fewer than 10 variants are not shown.

**Figure 5. ACMG/AMP criteria applied by CLAVAGEN and by expert panels.** (a) Precision and recall, with Wilson 95% confidence intervals, for criteria CLAVAGEN can evaluate (criteria observed at least 5 times). (b) Frequency of criteria applied by experts that require clinical, functional, segregation or case data.

**Figure 6. Causes of disagreement.** (a) Main cause per discordant variant by expert tier: missing evidence (the expert class recalculated from computable criteria equals CLAVAGEN's), criterion assignment, or evidence strength. (b) Most frequent expert-to-CLAVAGEN class changes.
"""


def make_figures(summary, results, out_dir):
    """Write all figures; returns the list of files. Raises ImportError if matplotlib is missing."""
    plt = _style()
    os.makedirs(out_dir, exist_ok=True)
    made = []
    for fn in (lambda: fig_selection(plt, summary, out_dir, made), lambda: fig_confusion(plt, summary, out_dir, made),
               lambda: fig_agreement(plt, summary, out_dir, made), lambda: fig_panels(plt, summary, out_dir, made),
               lambda: fig_criteria(plt, summary, out_dir, made), lambda: fig_disagreement(plt, results, out_dir, made)):
        fn()
    f = (summary.get("run") or {}).get("erepo_filter") or {}
    with open(os.path.join(out_dir, "figure_legends.md"), "w") as fh:
        fh.write(LEGENDS.format(total=f.get("total", 0), analysed=summary["analysed"]))
    made.append(os.path.join(out_dir, "figure_legends.md"))
    return made
