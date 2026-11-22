"""Pairwise alignment: affine-gap Smith-Waterman (local) and Needleman-Wunsch/Gotoh (global).

A pure-Python implementation is always available. If Biopython is installed it
is used automatically for speed (``backend='auto'``).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .genetic_code import blosum62
from .seqio import IUPAC_DNA

NEG = -1e18


@dataclass
class ScoringScheme:
    mol_type: str = "dna"
    match: float = 2.0
    mismatch: float = -3.0
    gap_open: float = 5.0     # penalty for a gap of length 1
    gap_extend: float = 2.0   # penalty per additional gap position

    @classmethod
    def default(cls, mol_type):
        if mol_type == "protein":
            return cls(mol_type="protein", gap_open=11.0, gap_extend=1.0)
        return cls(mol_type="dna")

    def score(self, a: str, b: str) -> float:
        if self.mol_type == "protein":
            return float(blosum62(a, b))
        if a == b and a in "ACGT":
            return self.match
        if a == "N" or b == "N":
            return 0.0
        sa, sb = IUPAC_DNA.get(a, "ACGT"), IUPAC_DNA.get(b, "ACGT")
        if set(sa) & set(sb):
            # IUPAC ambiguity (e.g. heterozygous call) compatible with reference
            return self.match * 0.75
        return self.mismatch


@dataclass
class Alignment:
    ref_id: str
    query_id: str
    mode: str
    aligned_ref: str
    aligned_query: str
    ref_start: int  # 0-based inclusive
    ref_end: int    # exclusive
    query_start: int
    query_end: int
    score: float
    ref_len: int
    query_len: int
    mol_type: str
    backend: str = "python"
    stats: dict = field(default_factory=dict)

    def compute_stats(self, scheme: ScoringScheme):
        L = len(self.aligned_ref)
        ident = sim = gaps = amb = 0
        for a, b in zip(self.aligned_ref, self.aligned_query):
            if a == "-" or b == "-":
                gaps += 1
            elif a == b:
                ident += 1
                sim += 1
            else:
                if scheme.mol_type == "dna" and b in IUPAC_DNA and a in IUPAC_DNA.get(b, ""):
                    amb += 1
                    sim += 1
                elif scheme.mol_type == "protein" and blosum62(a, b) > 0:
                    sim += 1
        self.stats = {
            "alignment_length": L,
            "identities": ident,
            "identity_pct": round(100.0 * ident / L, 2) if L else 0.0,
            "similarity_pct": round(100.0 * sim / L, 2) if L else 0.0,
            "gaps": gaps,
            "gap_pct": round(100.0 * gaps / L, 2) if L else 0.0,
            "ambiguous_compatible": amb,
            "ref_coverage_pct": round(100.0 * (self.ref_end - self.ref_start) / self.ref_len, 2) if self.ref_len else 0.0,
            "query_coverage_pct": round(100.0 * (self.query_end - self.query_start) / self.query_len, 2) if self.query_len else 0.0,
        }
        return self.stats

    # ------------------------------------------------------------------ text
    def match_line(self, scheme: ScoringScheme) -> str:
        out = []
        for a, b in zip(self.aligned_ref, self.aligned_query):
            if a == "-" or b == "-":
                out.append(" ")
            elif a == b:
                out.append("|")
            elif scheme.mol_type == "dna":
                out.append(":" if (b in IUPAC_DNA and a in IUPAC_DNA.get(b, "")) else ".")
            else:
                out.append("+" if blosum62(a, b) > 0 else ".")
        return "".join(out)

    def format_text(self, scheme: ScoringScheme, width: int = 60) -> str:
        ml = self.match_line(scheme)
        lines = []
        st = self.stats or self.compute_stats(scheme)
        lines.append(f"# Local alignment ({self.mode}, backend={self.backend})" if self.mode == "local"
                     else f"# Alignment ({self.mode}, backend={self.backend})")
        lines.append(f"# Reference: {self.ref_id} (length {self.ref_len})")
        lines.append(f"# Query:     {self.query_id} (length {self.query_len})")
        lines.append(f"# Score: {self.score:.1f}   Length: {st['alignment_length']}   "
                     f"Identity: {st['identities']}/{st['alignment_length']} ({st['identity_pct']}%)   "
                     f"Similarity: {st['similarity_pct']}%   Gaps: {st['gaps']} ({st['gap_pct']}%)")
        lines.append(f"# Reference span: {self.ref_start + 1}-{self.ref_end}   "
                     f"Query span: {self.query_start + 1}-{self.query_end}")
        lines.append("# Match line: '|' identical, ':' IUPAC-compatible (heterozygous), "
                     "'+' similar (protein), '.' mismatch, ' ' gap")
        lines.append("")
        rp, qp = self.ref_start, self.query_start
        nw = max(len(str(self.ref_end)), len(str(self.query_end)))
        for i in range(0, len(self.aligned_ref), width):
            rs = self.aligned_ref[i:i + width]
            qs = self.aligned_query[i:i + width]
            r_from, q_from = rp + 1, qp + 1
            rp += sum(1 for c in rs if c != "-")
            qp += sum(1 for c in qs if c != "-")
            lines.append(f"Ref   {r_from:>{nw}}  {rs}  {rp}")
            lines.append(f"      {'':>{nw}}  {ml[i:i + width]}")
            lines.append(f"Query {q_from:>{nw}}  {qs}  {qp}")
            lines.append("")
        return "\n".join(lines)


# --------------------------------------------------------------------- core

def _gotoh(ref: str, qry: str, scheme: ScoringScheme, local: bool):
    n, m = len(ref), len(qry)
    go, ge = scheme.gap_open, scheme.gap_extend
    prof = {c: [scheme.score(c, q) for q in qry] for c in set(ref)}

    tb = []
    if local:
        Mp = [0.0] * (m + 1)
        Xp = [NEG] * (m + 1)
        Yp = [NEG] * (m + 1)
        tb.append(bytearray(m + 1))
    else:
        Mp = [0.0] + [NEG] * m
        Xp = [NEG] * (m + 1)
        Yp = [NEG] + [-(go + (j - 1) * ge) for j in range(1, m + 1)]
        row0 = bytearray(m + 1)
        for j in range(2, m + 1):
            row0[j] = 2 << 4
        tb.append(row0)

    best, bi, bj = 0.0, 0, 0
    for i in range(1, n + 1):
        Mc = [NEG] * (m + 1)
        Xc = [NEG] * (m + 1)
        Yc = [NEG] * (m + 1)
        row = bytearray(m + 1)
        if not local:
            Xc[0] = -(go + (i - 1) * ge)
            row[0] = (0 if i == 1 else 1) << 2
        p = prof[ref[i - 1]]
        for j in range(1, m + 1):
            a, b, c = Mp[j - 1], Xp[j - 1], Yp[j - 1]
            if a >= b and a >= c:
                bm, mp = a, 0
            elif b >= c:
                bm, mp = b, 1
            else:
                bm, mp = c, 2
            if local and bm <= 0.0:
                bm, mp = 0.0, 3
            mv = bm + p[j - 1]
            # X: gap in query (consumes reference)
            a, b, c = Mp[j] - go, Xp[j] - ge, Yp[j] - go
            if a >= b and a >= c:
                xv, xp = a, 0
            elif b >= c:
                xv, xp = b, 1
            else:
                xv, xp = c, 2
            # Y: gap in reference (consumes query)
            a, b, c = Mc[j - 1] - go, Xc[j - 1] - go, Yc[j - 1] - ge
            if a >= b and a >= c:
                yv, yp = a, 0
            elif b >= c:
                yv, yp = b, 1
            else:
                yv, yp = c, 2
            Mc[j], Xc[j], Yc[j] = mv, xv, yv
            row[j] = mp | (xp << 2) | (yp << 4)
            if local and mv > best:
                best, bi, bj = mv, i, j
        tb.append(row)
        Mp, Xp, Yp = Mc, Xc, Yc

    if local:
        i, j, state, score = bi, bj, 0, best
    else:
        i, j = n, m
        vals = (Mp[m], Xp[m], Yp[m])
        state = max(range(3), key=lambda k: vals[k])
        score = vals[state]

    end_i, end_j = i, j
    ar, aq = [], []
    while i > 0 or j > 0:
        code = tb[i][j]
        if state == 0:
            if i == 0 or j == 0:
                break
            ptr = code & 3
            ar.append(ref[i - 1])
            aq.append(qry[j - 1])
            i -= 1
            j -= 1
            if ptr == 3:
                break
            state = ptr
        elif state == 1:
            ar.append(ref[i - 1])
            aq.append("-")
            state = (code >> 2) & 3
            i -= 1
        else:
            ar.append("-")
            aq.append(qry[j - 1])
            state = (code >> 4) & 3
            j -= 1
        if local and (i == 0 or j == 0) and state == 0:
            break
    return "".join(reversed(ar)), "".join(reversed(aq)), i, end_i, j, end_j, score


def _biopython_align(ref, qry, scheme, local):
    from Bio import Align  # noqa
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.mode = "local" if local else "global"
    if scheme.mol_type == "protein":
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    else:
        alpha = "".join(sorted(IUPAC_DNA))
        mat = substitution_matrices.Array(alphabet=alpha, dims=2)
        for x in alpha:
            for y in alpha:
                mat[x, y] = scheme.score(x, y)
        aligner.substitution_matrix = mat
    aligner.open_gap_score = -scheme.gap_open
    aligner.extend_gap_score = -scheme.gap_extend
    aln = aligner.align(ref, qry)[0]
    rb, qb = aln.aligned
    ar, aq = [], []
    for k in range(len(rb)):
        (rs, re_), (qs, qe) = rb[k], qb[k]
        if k > 0:
            prs, pqs = rb[k - 1][1], qb[k - 1][1]
            ar.append(ref[prs:rs]); aq.append("-" * (rs - prs))
            ar.append("-" * (qs - pqs)); aq.append(qry[pqs:qs])
        ar.append(ref[rs:re_]); aq.append(qry[qs:qe])
    r0, r1 = int(rb[0][0]), int(rb[-1][1])
    q0, q1 = int(qb[0][0]), int(qb[-1][1])
    sr, sq = "".join(ar), "".join(aq)
    if not local:
        # add terminal gaps
        sr = ref[:r0] + "-" * q0 + sr + ref[r1:] + "-" * (len(qry) - q1)
        sq = "-" * r0 + qry[:q0] + sq + "-" * (len(ref) - r1) + qry[q1:]
        r0, r1, q0, q1 = 0, len(ref), 0, len(qry)
    return sr, sq, r0, r1, q0, q1, float(aln.score)


def align(ref_id, ref, query_id, qry, mol_type="dna", mode="local", scheme=None, backend="auto",
          extend_ends=True, max_overhang=15):
    scheme = scheme or ScoringScheme.default(mol_type)
    local = mode == "local"
    used = "python"
    res = None
    if backend in ("auto", "biopython"):
        try:
            res = _biopython_align(ref, qry, scheme, local)
            used = "biopython"
        except Exception:
            if backend == "biopython":
                raise
            res = None
    if res is None:
        res = _gotoh(ref, qry, scheme, local)
        ar, aq, r0, r1, q0, q1, sc = res
    else:
        ar, aq, r0, r1, q0, q1, sc = res
    if local and extend_ends:
        ar, aq, r0, r1, q0, q1 = _extend_ends(ref, qry, ar, aq, r0, r1, q0, q1, scheme, max_overhang)
    aln = Alignment(ref_id, query_id, mode, ar, aq, r0, r1, q0, q1, sc,
                    len(ref), len(qry), mol_type, backend=used)
    aln.compute_stats(scheme)
    return aln


def _extend_ends(ref, qry, ar, aq, r0, r1, q0, q1, scheme, max_overhang=15):
    """Recover terminal query bases soft-clipped by local alignment (ungapped extension).

    Local alignment drops mismatches at the very ends of the query (e.g. a start-codon
    variant). Short overhangs (<= max_overhang) are re-attached ungapped if at least
    half of the extended positions match, giving semi-global (query-complete) behaviour.
    """
    def ok(rs, qs):
        if not rs:
            return False
        m = sum(1 for a, b in zip(rs, qs) if scheme.score(a, b) > 0)
        return m * 2 >= len(rs)
    k = q0
    if 0 < k <= max_overhang and r0 >= k and ok(ref[r0 - k:r0], qry[:k]):
        ar, aq = ref[r0 - k:r0] + ar, qry[:k] + aq
        r0, q0 = r0 - k, 0
    k = len(qry) - q1
    if 0 < k <= max_overhang and len(ref) - r1 >= k and ok(ref[r1:r1 + k], qry[q1:]):
        ar, aq = ar + ref[r1:r1 + k], aq + qry[q1:]
        r1, q1 = r1 + k, len(qry)
    return ar, aq, r0, r1, q0, q1
