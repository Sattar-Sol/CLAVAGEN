"""Extract variants from an alignment and normalise them (HGVS 3' rule)."""

from __future__ import annotations

from dataclasses import dataclass, field

from .seqio import HET_CODES, IUPAC_DNA


@dataclass
class RawVariant:
    ref_start: int          # 0-based inclusive on reference
    ref_end: int            # 0-based exclusive (== ref_start for insertions)
    ref: str                # reference allele ('' for insertion)
    alt: str                # alternate allele ('' for deletion)
    kind: str               # SNV | MNV | deletion | insertion | duplication | delins
    query_pos: int          # 0-based position in query
    zygosity: str = "homozygous/hemizygous (apparent)"
    repeat_unit: str = ""
    repeat_copies: int = 0
    notes: list = field(default_factory=list)

    @property
    def length_change(self):
        return len(self.alt) - len(self.ref)


def _min_period(s: str) -> str:
    for p in range(1, len(s) + 1):
        if len(s) % p == 0 and s[:p] * (len(s) // p) == s:
            return s[:p]
    return s


def _count_tract(ref: str, pos: int, unit: str) -> int:
    """Count tandem copies of ``unit`` in ref around ``pos`` (pos = left edge)."""
    if not unit:
        return 0
    k = len(unit)
    # move left to start of tract
    s = pos
    while s - k >= 0 and ref[s - k:s] == unit:
        s -= k
    n, e = 0, s
    while ref[e:e + k] == unit:
        n += 1
        e += k
    return n


def normalize(v: RawVariant, ref: str) -> RawVariant:
    """Shift pure indels to the most 3' position and detect duplications."""
    if v.kind == "deletion":
        s, e = v.ref_start, v.ref_end
        while e < len(ref) and ref[s] == ref[e]:
            s += 1
            e += 1
        v.ref_start, v.ref_end, v.ref = s, e, ref[s:e]
        unit = _min_period(v.ref)
        v.repeat_unit, v.repeat_copies = unit, _count_tract(ref, s, unit)
    elif v.kind in ("insertion", "duplication"):
        p, ins = v.ref_start, v.alt
        while p < len(ref) and ins[0] == ref[p]:
            ins = ins[1:] + ref[p]
            p += 1
        v.ref_start = v.ref_end = p
        v.alt = ins
        k = len(ins)
        if p - k >= 0 and ref[p - k:p] == ins:
            v.kind = "duplication"
        unit = _min_period(ins)
        v.repeat_unit, v.repeat_copies = unit, _count_tract(ref, p, unit) if p - len(unit) >= 0 else 0
    return v


def call_variants(aln, ref: str, mol_type: str = "dna", merge_adjacent: bool = True):
    """Return (variants, qc) from an Alignment object."""
    variants = []
    qc = {"no_call_positions": [], "ambiguous_uncallable": []}
    rp, qp = aln.ref_start, aln.query_start
    grp = None  # dict(start, ref, alt, has_s, has_i, has_d, qpos)

    def flush():
        nonlocal grp
        if not grp:
            return
        r, a = grp["ref"], grp["alt"]
        if grp["has_i"] and not grp["has_d"] and not grp["has_s"]:
            kind = "insertion"
        elif grp["has_d"] and not grp["has_i"] and not grp["has_s"]:
            kind = "deletion"
        elif not grp["has_i"] and not grp["has_d"]:
            kind = "SNV" if len(r) == 1 else "MNV"
        else:
            kind = "delins"
        if kind == "MNV" and mol_type == "protein" and not merge_adjacent:
            for k in range(len(r)):
                variants.append(RawVariant(grp["start"] + k, grp["start"] + k + 1, r[k], a[k], "SNV", grp["qpos"] + k))
        else:
            end = grp["start"] + len(r)
            variants.append(RawVariant(grp["start"], end, r, a, kind, grp["qpos"]))
        grp = None

    for rc, qc_ in zip(aln.aligned_ref, aln.aligned_query):
        if rc != "-" and qc_ != "-":
            if rc == qc_:
                flush()
            elif mol_type == "dna" and qc_ not in "ACGT":
                flush()
                if qc_ == "N":
                    qc["no_call_positions"].append(rp + 1)
                elif qc_ in HET_CODES and rc in HET_CODES[qc_]:
                    alt = HET_CODES[qc_].replace(rc, "")
                    variants.append(RawVariant(rp, rp + 1, rc, alt, "SNV", qp, zygosity="heterozygous",
                                               notes=[f"IUPAC code {qc_} in query interpreted as heterozygous {rc}/{alt}"]))
                else:
                    qc["ambiguous_uncallable"].append(f"{rp + 1}:{rc}>{qc_}({IUPAC_DNA.get(qc_, '?')})")
            elif mol_type == "protein" and qc_ == "X":
                flush()
                qc["no_call_positions"].append(rp + 1)
            else:
                if grp is None:
                    grp = dict(start=rp, ref="", alt="", has_s=False, has_i=False, has_d=False, qpos=qp)
                grp["ref"] += rc
                grp["alt"] += qc_
                grp["has_s"] = True
                if not merge_adjacent:
                    flush()
            rp += 1
            qp += 1
        elif rc == "-":
            if grp is None:
                grp = dict(start=rp, ref="", alt="", has_s=False, has_i=False, has_d=False, qpos=qp)
            grp["alt"] += qc_
            grp["has_i"] = True
            qp += 1
        else:
            if grp is None:
                grp = dict(start=rp, ref="", alt="", has_s=False, has_i=False, has_d=False, qpos=qp)
            grp["ref"] += rc
            grp["has_d"] = True
            rp += 1
    flush()
    variants = [normalize(v, ref) for v in variants]
    return variants, qc
