"""Truncating changes in protein-level queries (V11).

A frameshift changes every amino acid after the break and then stops early. A local alignment
of such a protein matches the reference only up to the break, so the changed tail and the early
stop fall outside the alignment and no variant is called. This module looks at what the
alignment left out and calls:

* a frameshift   - the query continues with residues unrelated to the reference (then stops);
* a nonsense     - the query stops ('*') right where the alignment ends;
* a stop-loss    - the query runs on past the end of the reference protein;

and, when a coding sequence is available, finds the small DNA insertion/deletion that explains
the observed frameshift tail.
"""
from .caller import RawVariant, normalize
from .consequence import hgvs_raw
from .genetic_code import translate

MIN_ALIGNED = 10        # aligned residues needed before a truncation is called
SIMILAR_TAIL = 0.5      # identity of the tail to the reference above which it is not a frameshift


def _identity(a, b):
    n = min(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x == y) / n if n else 0.0


def _last_solid_end(aln, run=4, max_trailing=8):
    """End (ref, query) of the last run of >= ``run`` identical residues, if only a short stretch
    containing a mismatch or gap follows it at the end of the alignment; otherwise None."""
    rp, qp = aln.ref_start, aln.query_start
    cur, last, last_end, imperfect_after = 0, None, None, False
    cols = list(zip(aln.aligned_ref, aln.aligned_query))
    for k, (a, b) in enumerate(cols):
        if a != "-":
            rp += 1
        if b != "-":
            qp += 1
        if a == b and a != "-":
            cur += 1
            if cur >= run:
                last, last_end, imperfect_after = (rp, qp), k, False
        else:
            cur = 0
            if last is not None:
                imperfect_after = True
    if last is None or not imperfect_after or len(cols) - 1 - last_end > max_trailing:
        return None
    return last


def detect_protein_truncation(aln, ref_prot, qseq, align_fn=None):
    """Return (RawVariant | None, notes, warnings) for the part of a protein query outside the alignment.

    When a frameshift is called, the RawVariant's ``ref_start`` may lie inside the aligned region
    (the first new residues can align by chance); callers should drop variants called from
    ``ref_start`` onwards.
    """
    notes, warns = [], []
    ref_prot = ref_prot.rstrip("*")
    q = qseq.upper()
    stop = q.find("*")
    core = q[:stop] if stop >= 0 else q
    has_stop = stop >= 0
    rend, qend, rlen = aln.ref_end, aln.query_end, len(ref_prot)
    if aln.ref_end - aln.ref_start < MIN_ALIGNED or qend > len(core) + 1:
        return None, notes, warns
    qend = min(qend, len(core))
    tail = core[qend:]

    # (d) query runs past the reference end: stop-loss / C-terminal extension
    if rend >= rlen:
        if tail:
            alt = tail + ("*" if has_stop else "")
            notes.append(f"Protein query continues {len(tail)} residue(s) beyond the reference stop: called as a "
                         "stop-loss (C-terminal extension). Exclude a tag or vector sequence.")
            return RawVariant(rlen, rlen, "", alt, "extension", qend), notes, warns
        return None, notes, warns

    # (a) nonsense directly after the aligned part
    if not tail and has_stop:
        notes.append(f"Protein query stops at residue {rend + 1} (reference length {rlen}): called as a premature stop.")
        return RawVariant(rend, rend + 1, ref_prot[rend], "*", "SNV", qend), notes, warns

    # (c) query simply ends: partial sequence or unmarked truncation
    if not tail:
        warns.append(f"Protein query ends at residue {rend} of {rlen}; residues {rend + 1}-{rlen} were not analysed. "
                     "If this is the complete protein it is truncated - end a complete sequence with '*' so a "
                     "premature stop can be called.")
        return None, notes, warns

    # (b) divergent tail: frameshift unless the tail is really the continuing reference.
    # The first new residues may have been aligned (as mismatches or chance matches): move the
    # break back to the end of the last solid run of identical residues.
    solid = _last_solid_end(aln)
    if solid is not None:
        rend, qend = solid
        tail = core[qend:]
    ref_next = ref_prot[rend:rend + len(tail)]
    if len(tail) >= 3 and _identity(tail, ref_next) >= SIMILAR_TAIL:
        warns.append(f"Residues {rend + 1}-{rend + len(tail)} differ partly from the reference but were not aligned; "
                     "check sequence quality in this region.")
        return None, notes, warns
    if align_fn is not None and len(tail) >= 8:
        sub = align_fn(ref_prot[rend:], tail)
        span = sub.query_end - sub.query_start
        if span >= 0.7 * len(tail) and sub.stats.get("identity_pct", 0) >= 70:
            warns.append(f"The query tail after residue {rend} matches the reference further downstream; the "
                         "alignment was split (possible large in-frame indel). Check the alignment.")
            return None, notes, warns
    alt = tail + ("*" if has_stop else "")
    stop_txt = (f"then stops after {len(tail)} new residue(s)" if has_stop
                else f"and {len(tail)} new residue(s) follow without a stop in the supplied sequence")
    notes.append(f"Protein query diverges from the reference after residue {rend} {stop_txt}: called as a probable "
                 "frameshift. The DNA change must be confirmed.")
    return RawVariant(rend, rlen, ref_prot[rend:], alt, "frameshift", qend), notes, warns


FS_DELETIONS = (1, 2, 4, 5)
FS_INSERTIONS = ("A", "C", "G", "T") + tuple(a + b for a in "ACGT" for b in "ACGT")


def infer_frameshift_dna(cds, utr3, tprot_off, v):
    """Find the frameshifting insertions/deletions in ``cds`` that give the observed protein.

    ``cds`` starts at c.1 (A of ATG); ``tprot_off`` is the offset of the reference protein in the
    transcript translation; ``v`` is the annotated variant (needs aa_pos and alt_tail).
    Returns a list of c. names (3'-normalised, de-duplicated).
    """
    i = tprot_off + v["aa_pos"] - 1                 # first changed codon (0-based, transcript)
    tail = v.get("alt_tail") or ""
    has_stop = tail.endswith("*")
    full = cds + (utr3 or "")
    expected_prefix = translate(cds[:i * 3])        # unchanged part
    lo, hi = max(0, (i - 3) * 3), min(len(cds), (i + 1) * 3 + 1)
    c0 = max(0, i - 3) * 3
    seen, names = set(), []
    edits = [("del", k, "") for k in FS_DELETIONS] + [("ins", len(s), s) for s in FS_INSERTIONS]
    for pos in range(lo, hi):
        for typ, k, ins in edits:
            if typ == "del":
                if pos + k > len(cds):
                    continue
                mut = full[:pos] + full[pos + k:]
            else:
                mut = full[:pos] + ins + full[pos:]
            prot = translate(mut[c0:], to_stop=True)
            got = translate(mut[:c0]) + prot
            if not got.startswith(expected_prefix):
                continue
            new = got[i:]
            if has_stop:
                if not new.endswith("*") or new != tail:
                    continue
            elif not new.startswith(tail):
                continue
            if mut in seen:
                continue
            seen.add(mut)
            rv = (RawVariant(pos, pos + k, full[pos:pos + k], "", "deletion", 0) if typ == "del"
                  else RawVariant(pos, pos, "", ins, "insertion", 0))
            rv = normalize(rv, full)
            name = hgvs_raw(rv, prefix="c.")
            if rv.kind == "insertion":
                name = f"c.{rv.ref_start}_{rv.ref_start + 1}ins{rv.alt}"
            if name not in names:
                names.append(name)
    return names
