"""SVG variant maps (no external dependencies)."""

from __future__ import annotations

from html import escape

from .consequence import CONSEQUENCE_LABEL

ACMG_COLORS = {
    "Pathogenic": "#b2182b", "Likely pathogenic": "#ef6548", "Uncertain significance": "#e3a008",
    "Likely benign": "#6baed6", "Benign": "#2166ac",
}
CONS_GROUPS = [
    ("Null / LOF", "#7a0177", {"stop_gained", "frameshift_variant", "splice_donor_variant",
                               "splice_acceptor_variant", "start_lost", "exon_loss_variant"}),
    ("Missense", "#e66101", {"missense_variant"}),
    ("In-frame / stop-loss", "#5e3c99", {"inframe_insertion", "inframe_deletion",
                                         "protein_altering_variant", "stop_lost"}),
    ("Splice region", "#1b9e77", {"splice_region_variant"}),
    ("Synonymous", "#4daf4a", {"synonymous_variant", "stop_retained_variant"}),
    ("UTR / intronic / other", "#7f7f7f", set()),
]
DOMAIN_COLORS = ["#a6cee3", "#b2df8a", "#fdbf6f", "#cab2d6", "#fb9a99", "#ffff99"]
FONT = "font-family:Helvetica,Arial,sans-serif"
MONO = "font-family:'DejaVu Sans Mono',Menlo,Consolas,monospace"


def cons_color(cons):
    for _, col, s in CONS_GROUPS:
        if cons in s:
            return col
    return CONS_GROUPS[-1][1]


def _t(x, y, s, size=11, anchor="start", color="#222", weight="normal", style=FONT, extra=""):
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" fill="{color}" '
            f'font-weight="{weight}" style="{style}" {extra}>{escape(str(s))}</text>')


def variant_map_svg(sample, width=1100):
    """Build the per-sequence variant map."""
    ref = sample["reference"]
    aln = sample["alignment"]
    tm = sample.get("tm")
    ctx = sample["context"]
    variants = sample["variants"]
    protein_mode = sample["mode"] == "protein"
    L = len(sample["ref_seq"])
    ml, mr = 70, 40
    span = width - ml - mr

    def X(p):  # p: 0-based coordinate (may be fractional)
        return ml + span * (p / max(L, 1))

    lane_h = 36
    top = 78
    # pre-assign label lanes so the figure height fits the number of lanes used
    max_lanes = 6
    lane_end = [-1e9] * max_lanes
    lane_of = {}
    for v in sorted(variants, key=lambda q: q["ref_pos"]):
        x = X(v["ref_pos"] - 0.5)
        lab1 = v["hgvs_c_short"] or v["hgvs_p"]
        lab2 = (v["hgvs_p"] + "  " if v["hgvs_c_short"] and v["hgvs_p"] else "") + CONSEQUENCE_LABEL.get(v["consequence"], v["consequence"])
        wlab = 6.3 * max(len(lab1) + 6, len(lab2)) + 16
        lane = next((k for k in range(max_lanes) if lane_end[k] < x - 4), min(range(max_lanes), key=lambda k: lane_end[k]))
        lane_end[lane] = x + wlab
        lane_of[v["id"]] = (lane, lab1, lab2)
    lanes = max([l for l, _, _ in lane_of.values()], default=0) + 1
    bar_y = top + lanes * lane_h + 24
    parts = []
    # ------------------------------------------------ header
    title = f"Variant map - {sample['query'].id} vs {ref.id}"
    sub = (f"{'Protein' if protein_mode else 'DNA'} | reference length {L} {'aa' if protein_mode else 'bp'} | "
           f"gene {ctx.get('gene', '-')} {ctx.get('transcript', '')} | identity {aln.stats['identity_pct']}% | "
           f"{len(variants)} variant(s)")
    parts.append(_t(ml, 26, title, 16, weight="bold"))
    parts.append(_t(ml, 46, sub, 11, color="#555"))

    # ------------------------------------------------ reference track
    if not protein_mode and tm is not None:
        parts.append(f'<line x1="{X(0)}" y1="{bar_y + 11}" x2="{X(L)}" y2="{bar_y + 11}" stroke="#999" stroke-width="2"/>')
        for s, e in tm.exons:
            parts.append(f'<rect x="{X(s):.1f}" y="{bar_y + 3}" width="{max(X(e) - X(s), 1):.1f}" height="16" '
                         f'fill="#d9e2ec" stroke="#8899aa" stroke-width="0.6"/>')
        # CDS blocks
        for s, e in tm.exons:
            cs, ce = max(s, tm.cds_start), min(e, tm.cds_end)
            if ce > cs:
                parts.append(f'<rect x="{X(cs):.1f}" y="{bar_y}" width="{max(X(ce) - X(cs), 1):.1f}" height="22" '
                             f'fill="#3b6ea5" stroke="#2c5282" stroke-width="0.6"/>')
        if len(tm.exons) == 1 and tm.junctions_c:
            for jc in tm.junctions_c:
                if 0 < jc <= len(tm.coding):
                    g = tm.coding[jc - 1] + 1
                    parts.append(f'<line x1="{X(g):.1f}" y1="{bar_y - 4}" x2="{X(g):.1f}" y2="{bar_y + 26}" '
                                 f'stroke="#fff" stroke-width="2" stroke-dasharray="3,2"/>')
                    parts.append(_t(X(g), bar_y + 36, f"ex{tm.junctions_c.index(jc) + 1}|ex{tm.junctions_c.index(jc) + 2}",
                                    8, "middle", "#666"))
        parts.append(_t(ml - 8, bar_y + 15, "Ref", 11, "end", weight="bold"))
    else:
        parts.append(f'<rect x="{X(0)}" y="{bar_y}" width="{span}" height="22" rx="4" fill="#3b6ea5"/>')
        parts.append(_t(ml - 8, bar_y + 15, "Prot", 11, "end", weight="bold"))

    # ------------------------------------------------ domain track
    dom_y = bar_y + 44
    doms = ctx.get("domains") or []
    row_end, drow = [], []
    for d in sorted(doms, key=lambda d: int(d["start"])):
        r = next((k for k, e in enumerate(row_end) if e < int(d["start"])), None)
        if r is None:
            row_end.append(int(d["end"]))
            r = len(row_end) - 1
        else:
            row_end[r] = int(d["end"])
        drow.append((d, r))
    n_drows = max(1, len(row_end))
    for i, (d, r) in enumerate(drow):
        s, e = int(d["start"]), int(d["end"])
        if protein_mode:
            x0, x1 = X(s - 1), X(e)
        elif tm is not None and tm.coding:
            c0 = min(max((s - 1) * 3, 0), len(tm.coding) - 1)
            c1 = min(e * 3 - 1, len(tm.coding) - 1)
            x0, x1 = X(tm.coding[c0]), X(tm.coding[c1] + 1)
        else:
            continue
        col = DOMAIN_COLORS[i % len(DOMAIN_COLORS)]
        stroke = "#b2182b" if d.get("critical") else "#666"
        yy = dom_y + r * 17
        parts.append(f'<rect x="{x0:.1f}" y="{yy}" width="{max(x1 - x0, 2):.1f}" height="14" rx="3" fill="{col}" '
                     f'stroke="{stroke}" stroke-width="{1.5 if d.get("critical") else 0.6}"/>')
        parts.append(_t((x0 + x1) / 2, yy + 11, d["name"] + (" *" if d.get("critical") else ""), 9, "middle"))
    if doms:
        parts.append(_t(ml - 8, dom_y + 11, "Domains", 10, "end"))

    # ------------------------------------------------ query coverage track
    cov_y = dom_y + (n_drows * 17 + 9 if doms else 0)
    parts.append(_t(ml - 8, cov_y + 9, "Query", 10, "end"))
    parts.append(f'<rect x="{X(aln.ref_start):.1f}" y="{cov_y}" width="{max(X(aln.ref_end) - X(aln.ref_start), 1):.1f}" '
                 f'height="10" fill="#9ccfa0" stroke="#5a9a5e" stroke-width="0.6"/>')
    nocalls = sample["qc"].get("no_call_positions", [])
    for p in nocalls:
        parts.append(f'<rect x="{X(p - 1):.1f}" y="{cov_y}" width="{max(X(p) - X(p - 1), 1):.1f}" height="10" fill="#555"/>')

    # ------------------------------------------------ axis
    ax_y = cov_y + 26
    parts.append(f'<line x1="{X(0)}" y1="{ax_y}" x2="{X(L)}" y2="{ax_y}" stroke="#444"/>')
    step = _nice_step(L)
    p = 0
    while p <= L:
        parts.append(f'<line x1="{X(p):.1f}" y1="{ax_y}" x2="{X(p):.1f}" y2="{ax_y + 5}" stroke="#444"/>')
        parts.append(_t(X(p), ax_y + 17, max(p, 1), 9, "middle", "#444"))
        p += step
    parts.append(_t(X(L), ax_y + 30, "reference position (aa)" if protein_mode else "reference position (nt)",
                    9, "end", "#666"))

    # ------------------------------------------------ lollipops
    for v in sorted(variants, key=lambda q: q["ref_pos"]):
        x = X(v["ref_pos"] - 0.5)
        lane, lab1, lab2 = lane_of[v["id"]]
        hy = bar_y - 22 - lane * lane_h
        cc = cons_color(v["consequence"])
        ac = ACMG_COLORS.get(v["classification"], "#999")
        parts.append(f'<line x1="{x:.1f}" y1="{bar_y}" x2="{x:.1f}" y2="{hy}" stroke="{cc}" stroke-width="1.6"/>')
        parts.append(f'<polygon points="{x - 5:.1f},{bar_y + 30} {x + 5:.1f},{bar_y + 30} {x:.1f},{bar_y + 23}" fill="{cc}"/>')
        parts.append(f'<circle cx="{x:.1f}" cy="{hy}" r="8" fill="{ac}" stroke="{cc}" stroke-width="3">'
                     f'<title>{escape(v["id"])}: {escape(lab1)} {escape(v["hgvs_p"])} | {escape(v["classification"])} '
                     f'| {escape(", ".join(v["criteria_met"]))}</title></circle>')
        parts.append(_t(x, hy + 3.5, v["class_abbr"], 7, "middle", "#fff", "bold"))
        zyg = " (het)" if v["zygosity"] == "heterozygous" else ""
        parts.append(_t(x + 12, hy - 2, lab1 + zyg, 10.5, weight="bold"))
        parts.append(_t(x + 12, hy + 10, lab2, 9.5, color=cc))

    # ------------------------------------------------ legend
    lg_y = ax_y + 50
    parts.append(_t(ml, lg_y, "ACMG class (circle fill):", 10, weight="bold"))
    xx = ml + 150
    for k, c in ACMG_COLORS.items():
        parts.append(f'<circle cx="{xx}" cy="{lg_y - 4}" r="6" fill="{c}"/>')
        parts.append(_t(xx + 10, lg_y, k, 10))
        xx += 22 + 6.2 * len(k)
    lg_y += 20
    parts.append(_t(ml, lg_y, "Variant type (stem/ring):", 10, weight="bold"))
    xx = ml + 150
    for name, col, _ in CONS_GROUPS:
        parts.append(f'<rect x="{xx - 6}" y="{lg_y - 10}" width="12" height="10" fill="{col}"/>')
        parts.append(_t(xx + 10, lg_y, name, 10))
        xx += 24 + 6.0 * len(name)
    lg_y += 20
    extra = "Dark-blue = CDS, light = UTR/exon, line = intron; dashed white = exon junction; green = query coverage; grey = no-call (N); * = critical domain"
    parts.append(_t(ml, lg_y, extra, 9, color="#666"))

    # ------------------------------------------------ zoom panels
    y = lg_y + 26
    zoom = []
    for v in variants[:16]:
        zoom.append(_zoom_panel(v, sample, protein_mode, tm))
    col_w = (width - ml) / 2
    for i, (h, body) in enumerate(zoom):
        cx = ml - 40 + (i % 2) * col_w
        if i % 2 == 0 and i:
            y += prev_h + 12
        if i % 2 == 0:
            prev_h = h
        else:
            prev_h = max(prev_h, h)
        parts.append(f'<g transform="translate({cx:.1f},{y:.1f})">{body}</g>')
    height = (y + (prev_h if zoom else 0) + 20)
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{int(height)}" '
           f'viewBox="0 0 {width} {int(height)}"><rect width="100%" height="100%" fill="#fff"/>'
           + "".join(parts) + "</svg>")
    return svg


def _nice_step(L):
    for s in (1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 5000, 10000, 20000, 50000):
        if L / s <= 12:
            return s
    return 100000


def _zoom_panel(v, sample, protein_mode, tm, flank=12):
    seq = sample["ref_seq"]
    s, e = v["ref_pos"] - 1, v["ref_end"]
    if v["kind"] in ("insertion", "duplication"):
        e = s
    a, b = max(0, s - flank), min(len(seq), e + flank)
    ref_w = seq[a:b]
    alt_ins = v["alt_allele"]
    if v["zygosity"] == "heterozygous" and len(v["ref_allele"]) == 1:
        from .seqio import IUPAC_DNA
        code = {vv: k for k, vv in IUPAC_DNA.items()}.get("".join(sorted(v["ref_allele"] + v["alt_allele"])), alt_ins)
        alt_ins = code
    alt_w = seq[a:s] + alt_ins + seq[e:b]
    cw = 9.0
    x0 = 70
    parts = []
    cc = cons_color(v["consequence"])
    ac = ACMG_COLORS.get(v["classification"], "#999")
    parts.append(f'<rect x="0" y="0" width="500" height="96" rx="6" fill="#fafafa" stroke="#ddd"/>')
    parts.append(f'<rect x="0" y="0" width="6" height="96" rx="2" fill="{ac}"/>')
    head = f"{v['id']}  {v['hgvs_c_short'] or v['hgvs_p']}  {v['hgvs_p'] if v['hgvs_c_short'] else ''}"
    parts.append(_t(14, 16, head, 11, weight="bold"))
    parts.append(_t(486, 16, f"{v['classification']} ({v['points']} pts)", 10, "end", ac, "bold"))
    parts.append(_t(14, 30, f"{CONSEQUENCE_LABEL.get(v['consequence'], v['consequence'])} | {v['zygosity']} | "
                            f"{', '.join(v['criteria_met']) or 'no criteria met'}", 9, color=cc))
    # ref / alt rows aligned with gaps
    n_left = s - a
    ra = ref_w[:n_left] + ref_w[n_left:n_left + (e - s)] + "-" * max(0, len(alt_ins) - (e - s)) + ref_w[n_left + (e - s):]
    aa_ = alt_w[:n_left] + alt_ins + "-" * max(0, (e - s) - len(alt_ins)) + alt_w[n_left + len(alt_ins):]
    hl0, hl1 = n_left, n_left + max(e - s, len(alt_ins))
    max_cols = 40                       # keep long changes (frameshift tails, large indels) inside the card
    if max(len(ra), len(aa_)) > max_cols:
        ra = ra[:max_cols - 1] + "…"
        aa_ = aa_[:max_cols - 1] + "…"
        hl1 = min(hl1, max_cols - 1)
        b = a + sum(1 for ch in ra[:-1] if ch != "-")
    parts.append(_t(14, 52, "Ref", 10, weight="bold"))
    parts.append(_t(14, 70, "Query", 10, weight="bold"))
    for k, ch in enumerate(ra):
        xx = x0 + k * cw
        if hl0 <= k < hl1:
            parts.append(f'<rect x="{xx - 1:.1f}" y="40" width="{cw:.1f}" height="36" fill="{cc}" opacity="0.18"/>')
        parts.append(_t(xx + cw / 2 - 1, 52, ch, 11, "middle", style=MONO))
    for k, ch in enumerate(aa_):
        xx = x0 + k * cw
        col = cc if hl0 <= k < hl1 else "#222"
        parts.append(_t(xx + cw / 2 - 1, 70, ch, 11, "middle", col, "bold" if hl0 <= k < hl1 else "normal", MONO))
    pos_label = f"{'aa' if protein_mode else 'ref'} {a + 1}-{b}"
    parts.append(_t(x0 + len(ra) * cw + 6, 61, pos_label, 8, color="#888"))
    # codon / aa line
    if not protein_mode and tm is not None and v.get("codon_change"):
        parts.append(_t(14, 88, f"Codon {v['codon_change']}   {v.get('aa_change', '')}   "
                                f"BLOSUM62 {v.get('blosum62', '-')}  Grantham {v.get('grantham', '-')}", 9, color="#444"))
    elif protein_mode and v.get("aa_change"):
        parts.append(_t(14, 88, f"{v.get('aa_change')}   BLOSUM62 {v.get('blosum62', '-')}  "
                                f"Grantham {v.get('grantham', '-')}", 9, color="#444"))
    elif v.get("region"):
        parts.append(_t(14, 88, f"Region: {v['region']}" + (f"  exon {v['exon']}" if v.get("exon") else ""), 9, color="#444"))
    return 96, "".join(parts)
