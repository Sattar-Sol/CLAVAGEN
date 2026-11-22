"""Sanger trace (.ab1 / ABIF) support.

* ABIF reader (Applied Biosystems format, standard library only)
* quality trimming (modified Mott algorithm on the basecaller quality values)
* heterozygous base calling from secondary peaks -> IUPAC consensus
* heterozygous indel detection: after a heterozygous insertion/deletion the two alleles are
  out of register and almost every peak is double. The trace is split into the reference-
  matching allele and the second allele, which is then aligned to the reference.
* chromatogram drawing (SVG) around any reference position
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from html import escape

from .genetic_code import COMPLEMENT

PAIR = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M"}
CHANNEL_COLORS = {"A": "#1a9850", "C": "#2b6cb0", "G": "#222222", "T": "#d7301f"}


# ------------------------------------------------------------------ ABIF reader
def is_abif(data: bytes) -> bool:
    return data[:4] == b"ABIF"


def read_abif(data: bytes) -> dict:
    """Return {TAGNUM: value} for the tags of an ABIF file (e.g. 'PBAS2', 'DATA9')."""
    if not is_abif(data):
        raise ValueError("Not an ABIF (.ab1) file")
    # header directory entry starts at byte 6
    _, _, _, _, n_el, _, dir_off, _ = struct.unpack(">4sihhiiii", data[6:34])
    tags = {}
    for i in range(n_el):
        e = data[dir_off + 28 * i: dir_off + 28 * (i + 1)]
        name, num, etype, esize, count, dsize, doff, _ = struct.unpack(">4sihhiiii", e)
        raw = e[20:20 + dsize] if dsize <= 4 else data[doff:doff + dsize]
        key = name.decode("latin-1") + str(num)
        try:
            if etype == 2:            # char array
                val = raw.decode("latin-1")
            elif etype == 4:          # short
                val = list(struct.unpack(f">{count}h", raw[:2 * count]))
            elif etype == 5:          # long
                val = list(struct.unpack(f">{count}i", raw[:4 * count]))
            elif etype == 1:          # byte
                val = list(raw)
            elif etype == 18:         # pString
                val = raw[1:1 + raw[0]].decode("latin-1") if raw else ""
            elif etype == 19:         # cString
                val = raw.split(b"\0")[0].decode("latin-1")
            elif etype == 7:
                val = list(struct.unpack(f">{count}f", raw[:4 * count]))
            else:
                val = raw
        except struct.error:
            val = raw
        if isinstance(val, list) and len(val) == 1 and etype in (4, 5, 7):
            val = val[0]
        tags[key] = val
    return tags


@dataclass
class Trace:
    name: str
    order: str                 # channel order of DATA9..12, e.g. "GATC"
    channels: dict             # base -> list of intensities
    peaks: list                # peak positions (one per base call)
    calls: str                 # basecaller calls
    quals: list                # basecaller quality values
    primary: str = ""
    secondary: str = ""
    ratio: list = field(default_factory=list)
    mixed: list = field(default_factory=list)
    trim: tuple = (0, 0)       # kept peak range [start, end)
    reversed: bool = False
    het_indel_start: int = None  # index (in trimmed coordinates) where the two alleles go out of register
    het_indel_end: int = None    # set when the out-of-register block is at the START of the read (reverse read
                                 # turned round): double peaks run from het_indel_start up to this index
    notes: list = field(default_factory=list)

    # trimmed views ---------------------------------------------------------
    def consensus(self):
        a, b = self.trim
        out = []
        for i in range(a, b):
            p, s = self.primary[i], self.secondary[i]
            out.append(PAIR.get("".join(sorted(p + s)), p) if self.mixed[i] and s != p else p)
        return "".join(out)

    def reverse_complement(self):
        n = len(self.peaks)
        L = len(next(iter(self.channels.values())))
        comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
        self.channels = {comp[b]: list(reversed(v)) for b, v in self.channels.items()}
        self.peaks = [L - 1 - p for p in reversed(self.peaks)]
        self.calls = self.calls[::-1].translate(COMPLEMENT)
        self.quals = list(reversed(self.quals))
        self.primary = self.primary[::-1].translate(COMPLEMENT)
        self.secondary = self.secondary[::-1].translate(COMPLEMENT)
        self.ratio = list(reversed(self.ratio))
        self.mixed = list(reversed(self.mixed))
        a, b = self.trim
        self.trim = (n - b, n - a)
        n_trim = b - a
        s0, e0 = self.het_indel_start, self.het_indel_end
        # A heterozygous indel puts the two alleles out of register from the indel to the END of the read
        # (in sequencing direction). Turning the read round moves that block to the START; keep it.
        if s0 is not None and e0 is None:
            self.het_indel_start, self.het_indel_end = 0, n_trim - s0
            self.notes.append("Reverse read turned round: the double-peak block of the heterozygous indel now lies "
                              "before the indel.")
        elif s0 is not None:
            self.het_indel_start, self.het_indel_end = n_trim - e0, None
        else:
            self.het_indel_start, self.het_indel_end = None, None
        self.reversed = not self.reversed

    @property
    def het_indel_site(self):
        """Trimmed index of the indel itself (where the reads go in or out of register)."""
        if self.het_indel_start is None:
            return None
        return self.het_indel_end if self.het_indel_end is not None else self.het_indel_start


def load_trace(data: bytes, name="trace", het_ratio=0.30, min_qual_trim=0.02):
    t = read_abif(data)
    order = (t.get("FWO_1") or "GATC")[:4]
    chans = {}
    for i, b in enumerate(order):
        d = t.get(f"DATA{9 + i}") or t.get(f"DATA{1 + i}")
        if d is None:
            raise ValueError("Trace data channels (DATA9-12) not found in .ab1 file")
        chans[b] = d if isinstance(d, list) else [d]
    calls = t.get("PBAS2") or t.get("PBAS1") or ""
    peaks = t.get("PLOC2") or t.get("PLOC1") or []
    if isinstance(peaks, int):
        peaks = [peaks]
    q = t.get("PCON2") or t.get("PCON1") or []
    quals = [ord(c) for c in q] if isinstance(q, str) else list(q)
    if len(quals) < len(calls):
        quals = quals + [20] * (len(calls) - len(quals))
    n = min(len(calls), len(peaks))
    tr = Trace(name=t.get("SMPL1") or name, order=order, channels=chans, peaks=peaks[:n], calls=calls[:n].upper(),
               quals=quals[:n])
    call_bases(tr, het_ratio)
    if not tr.quals or max(tr.quals) == 0:
        import math
        tr.quals = [int(min(50, max(0, -10 * math.log10(max(r, 1e-3))))) for r in tr.ratio]
        tr.notes.append("No basecaller quality values in file - quality estimated from peak purity.")
    # genuine double peaks (heterozygous positions) get low basecaller quality but are good signal:
    # do not let them be trimmed away
    q_eff = [max(q, 30) if m else q for q, m in zip(tr.quals, tr.mixed)]
    tr.trim = mott_trim(q_eff, min_qual_trim)
    if tr.trim == (0, len(tr.quals)) and sum(tr.quals) / max(1, len(tr.quals)) < 15:
        tr.notes.append("Low-quality trace: no high-quality region found; whole read used - interpret with caution.")
    tr.het_indel_start = detect_het_indel(tr)
    return tr


def _height(ch, pos):
    lo, hi = max(0, pos - 1), min(len(ch), pos + 2)
    return max(ch[lo:hi]) if hi > lo else 0


def _is_peak(ch, pos, w=2):
    """True if the channel has a genuine local maximum at this position (not a neighbour's tail)."""
    lo, hi = max(0, pos - w), min(len(ch), pos + w + 1)
    seg = ch[lo:hi]
    if len(seg) < 3:
        return False
    j = seg.index(max(seg))
    return 0 < j < len(seg) - 1


def call_bases(tr: Trace, het_ratio=0.30):
    prim, sec, ratio, mixed = [], [], [], []
    for i, pk in enumerate(tr.peaks):
        h = sorted(((_height(tr.channels[b], pk), b) for b in "ACGT"), reverse=True)
        (h1, b1), (h2, b2) = h[0], h[1]
        r = (h2 / h1) if h1 > 0 else 0.0
        call = tr.calls[i] if i < len(tr.calls) and tr.calls[i] in "ACGT" else b1
        if call != b1 and r >= 0.85:
            b1, b2 = call, b1 if call != b1 else b2
        prim.append(b1)
        sec.append(b2)
        ratio.append(round(r, 3))
        mixed.append(r >= het_ratio and h1 > 50 and _is_peak(tr.channels[b2], pk))
    tr.primary, tr.secondary, tr.ratio, tr.mixed = "".join(prim), "".join(sec), ratio, mixed


def mott_trim(quals, cutoff=0.05):
    """Modified Mott trimming: maximum-scoring segment of (cutoff - error probability)."""
    best = cur = 0.0
    start = bs = be = 0
    for i, q in enumerate(quals):
        cur += cutoff - 10 ** (-q / 10.0)
        if cur < 0:
            cur, start = 0.0, i + 1
        elif cur > best:
            best, bs, be = cur, start, i + 1
    if be - bs < 20:
        return (0, len(quals))
    return (bs, be)


def detect_het_indel(tr: Trace, window=20, frac=0.6, min_len=25):
    """Start (trimmed index) of a sustained run of double peaks, typical of a heterozygous indel."""
    a, b = tr.trim
    m = tr.mixed[a:b]
    n = len(m)
    for i in range(0, max(0, n - min_len)):
        w = m[i:i + window]
        if sum(w) >= frac * len(w) and sum(m[i:i + min_len]) >= frac * min_len and m[i]:
            # require the run to continue for most of the remaining read
            rest = m[i:]
            if sum(rest) >= 0.5 * len(rest):
                return i
    return None


def deconvolute(tr: Trace, ref_base_at):
    """Split a trace with a heterozygous indel into two alleles.

    ``ref_base_at(i)`` gives the reference base aligned to trimmed position i (or None).
    Returns (allele1, allele2): allele1 follows the reference in the mixed region (IUPAC kept
    upstream), allele2 carries the other base at every double peak.
    """
    a, b = tr.trim
    cons = tr.consensus()
    s0 = tr.het_indel_start or 0
    e0 = tr.het_indel_end if tr.het_indel_end is not None else len(cons)
    a1, a2 = list(cons), list(cons)
    for i in range(len(cons)):
        k = a + i
        p, s = tr.primary[k], tr.secondary[k]
        if i < s0 or i >= e0:
            # outside the double-peak block: allele 2 equals allele 1 except at heterozygous SNVs
            # (kept as IUPAC in allele 1)
            a2[i] = p
            continue
        r = ref_base_at(i)
        if tr.mixed[k] and r in (p, s):
            a1[i] = r
            a2[i] = s if r == p else p
        else:
            a1[i] = p
            a2[i] = p
    return "".join(a1), "".join(a2)


# ------------------------------------------------------------------ drawing
def chromatogram_svg(tr: Trace, center_idx, ref_window="", title="", flank=12, width=520, height=170):
    """Chromatogram around trimmed index ``center_idx`` (peak index = trim start + center_idx)."""
    a, _ = tr.trim
    k = a + center_idx
    lo, hi = max(0, k - flank), min(len(tr.peaks) - 1, k + flank)
    x0, x1 = tr.peaks[lo] - 6, tr.peaks[hi] + 6
    L = len(next(iter(tr.channels.values())))
    x0, x1 = max(0, x0), min(L - 1, x1)
    ymax = max(1, max(max(tr.channels[c][x0:x1 + 1] or [1]) for c in "ACGT"))
    top, bottom = 34, height - 40
    sx = lambda x: 10 + (x - x0) * (width - 20) / max(1, x1 - x0)  # noqa: E731
    sy = lambda y: bottom - (bottom - top) * y / ymax  # noqa: E731
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="#fff"/>',
             f'<text x="10" y="14" font-size="11" font-weight="bold" font-family="Helvetica,Arial">{escape(title)}</text>',
             f'<text x="{width - 10}" y="14" font-size="9.5" text-anchor="end" fill="#777" font-family="Helvetica,Arial">'
             f'{escape(tr.name)}{" (reverse read, shown reverse-complemented)" if tr.reversed else ""}</text>']
    cx = sx(tr.peaks[k])
    parts.append(f'<rect x="{cx - 9}" y="{top - 4}" width="18" height="{bottom - top + 34}" fill="#f0b400" opacity="0.18"/>')
    for c in "ACGT":
        pts = " ".join(f"{sx(x):.1f},{sy(tr.channels[c][x]):.1f}" for x in range(x0, x1 + 1, 1))
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{CHANNEL_COLORS[c]}" stroke-width="1.2"/>')
    for i in range(lo, hi + 1):
        x = sx(tr.peaks[i])
        call = tr.primary[i]
        if tr.mixed[i] and tr.secondary[i] != call:
            call = PAIR.get("".join(sorted(call + tr.secondary[i])), call)
        col = CHANNEL_COLORS.get(call, "#b45309")
        parts.append(f'<text x="{x:.1f}" y="{top - 8}" font-size="10" text-anchor="middle" font-weight="bold" '
                     f'fill="{col}" font-family="Menlo,Consolas,monospace">{call}</text>')
        j = i - k + flank
        if ref_window and 0 <= j < len(ref_window):
            parts.append(f'<text x="{x:.1f}" y="{bottom + 14}" font-size="10" text-anchor="middle" fill="#555" '
                         f'font-family="Menlo,Consolas,monospace">{ref_window[j]}</text>')
    parts.append(f'<line x1="10" y1="{bottom}" x2="{width - 10}" y2="{bottom}" stroke="#ccc"/>')
    parts.append(f'<text x="10" y="{bottom + 28}" font-size="9" fill="#777" font-family="Helvetica,Arial">'
                 'top: trace calls (IUPAC = double peak) | below line: reference | A green, C blue, G black, T red</text>')
    parts.append("</svg>")
    return "".join(parts)


# ------------------------------------------------------------------ synthetic traces (tests / demos)
def write_synthetic_ab1(allele1, allele2=None, path=None, name="synthetic", spacing=12, seed=1):
    """Write a minimal valid ABIF file simulating a Sanger read of one or two alleles
    (allele2 defaults to allele1 = homozygous). Used for tests and demonstrations."""
    import math
    import random
    rnd = random.Random(seed)
    allele2 = allele2 if allele2 is not None else allele1
    n = max(len(allele1), len(allele2))
    L = n * spacing + 30
    chans = {b: [rnd.randint(2, 15) for _ in range(L)] for b in "GATC"}
    for allele in (allele1, allele2):
        for i, base in enumerate(allele):
            if base not in "ACGT":
                continue
            c = 15 + i * spacing
            h = 520 + rnd.randint(-60, 60)
            for x in range(max(0, c - 8), min(L, c + 9)):
                chans[base][x] += int(h * math.exp(-((x - c) ** 2) / 7.0))
    peaks = [15 + i * spacing for i in range(n)]
    calls, quals = [], []
    for pk in peaks:
        h = sorted(((chans[b][pk], b) for b in "ACGT"), reverse=True)
        calls.append(h[0][1])
        r = h[1][0] / max(1, h[0][0])
        quals.append(max(5, min(60, int(-10 * math.log10(max(r, 1e-6))))))
    entries = [("FWO_", 1, 2, 1, b"GATC")]
    for i, b in enumerate("GATC"):
        entries.append(("DATA", 9 + i, 4, 2, struct.pack(f">{L}h", *[min(32767, v) for v in chans[b]])))
    entries += [("PLOC", 2, 4, 2, struct.pack(f">{n}h", *peaks)), ("PBAS", 2, 2, 1, "".join(calls).encode()),
                ("PCON", 2, 2, 1, bytes(quals)), ("SMPL", 1, 18, 1, bytes([len(name)]) + name.encode())]
    blob, dirs = b"", []
    data_start = 128
    for tag, num, etype, esize, raw in entries:
        count = len(raw) // esize if etype != 18 else len(raw)
        if len(raw) <= 4:
            off_field = raw.ljust(4, b"\0")
            dirs.append(struct.pack(">4sihhii", tag.encode(), num, etype, esize, count, len(raw)) + off_field +
                        struct.pack(">i", 0))
        else:
            off = data_start + len(blob)
            blob += raw
            dirs.append(struct.pack(">4sihhiiii", tag.encode(), num, etype, esize, count, len(raw), off, 0))
    dir_off = data_start + len(blob)
    header = b"ABIF" + struct.pack(">h", 101) + struct.pack(">4sihhiiii", b"tdir", 1, 1023, 28, len(dirs),
                                                             28 * len(dirs), dir_off, 0)
    out = header.ljust(data_start, b"\0") + blob + b"".join(dirs)
    if path:
        with open(path, "wb") as fh:
            fh.write(out)
    return out
