"""Transcript model, HGVS nomenclature and molecular consequence prediction."""

from __future__ import annotations

import bisect

from .genetic_code import AA3, find_orfs, translate

LOF_CONSEQUENCES = {"stop_gained", "frameshift_variant", "splice_donor_variant",
                    "splice_acceptor_variant", "start_lost", "exon_loss_variant"}

CONSEQUENCE_SEVERITY = [
    "exon_loss_variant", "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
    "frameshift_variant", "stop_lost", "start_lost", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant", "missense_variant", "splice_region_variant",
    "stop_retained_variant", "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "intron_variant", "upstream_gene_variant", "downstream_gene_variant", "non_coding_variant",
]

CONSEQUENCE_LABEL = {
    "exon_loss_variant": "Exon/boundary loss", "splice_acceptor_variant": "Splice acceptor",
    "splice_donor_variant": "Splice donor", "stop_gained": "Nonsense",
    "frameshift_variant": "Frameshift", "stop_lost": "Stop-loss", "start_lost": "Start-loss",
    "inframe_insertion": "In-frame insertion", "inframe_deletion": "In-frame deletion",
    "protein_altering_variant": "In-frame delins", "missense_variant": "Missense",
    "splice_region_variant": "Splice region", "stop_retained_variant": "Stop retained",
    "synonymous_variant": "Synonymous", "5_prime_UTR_variant": "5' UTR",
    "3_prime_UTR_variant": "3' UTR", "intron_variant": "Intronic",
    "upstream_gene_variant": "Upstream", "downstream_gene_variant": "Downstream",
    "non_coding_variant": "Non-coding",
}


def a3(c):
    return AA3.get(c, "Xaa")


# ===================================================================== model
class TranscriptModel:
    """Maps reference coordinates to HGVS c. coordinates.

    ``exons``: list of (start0, end0_exclusive) on the reference (genomic-style
    reference). If omitted, the whole reference is one exon (cDNA/mRNA reference).
    ``cds``: (start0, end0_exclusive) on the reference, including the stop codon.
    ``exon_boundaries_c``: for cDNA references, c. positions of the last base of
    each exon (except the last exon) - used for NMD and splice-region logic.
    """

    def __init__(self, ref_seq, cds=None, exons=None, exon_boundaries_c=None):
        self.ref = ref_seq
        self.L = len(ref_seq)
        self.warnings = []
        self.exons = sorted(tuple(e) for e in exons) if exons else [(0, self.L)]
        self.exonic = [g for s, e in self.exons for g in range(s, e)]
        self.eidx = {g: i for i, g in enumerate(self.exonic)}
        if cds is None:
            cds = self._auto_cds()
        self.cds_start, self.cds_end = cds
        self.coding = [g for s, e in self.exons for g in range(max(s, self.cds_start), min(e, self.cds_end))]
        self.cidx = {g: i for i, g in enumerate(self.coding)}
        self.cds_seq = "".join(ref_seq[g] for g in self.coding)
        if not self.coding:
            self.warnings.append("No coding sequence could be defined; consequences limited to positional terms.")
            self.cds_e0 = self.cds_e1 = None
        else:
            self.cds_e0 = self.eidx[self.coding[0]]
            self.cds_e1 = self.eidx[self.coding[-1]]
        self.utr3_seq = "".join(ref_seq[g] for g in self.exonic[self.cds_e1 + 1:]) if self.coding else ""
        self.protein = translate(self.cds_seq, to_stop=False)
        if self.coding:
            if len(self.cds_seq) % 3:
                self.warnings.append(f"CDS length {len(self.cds_seq)} is not a multiple of 3.")
            if not self.cds_seq.startswith("ATG"):
                self.warnings.append("CDS does not start with ATG.")
            if "*" in self.protein[:-1]:
                self.warnings.append("Internal stop codon(s) in reference CDS - check CDS coordinates.")
            if not self.protein.endswith("*"):
                self.warnings.append("Reference CDS does not end with a stop codon.")
        # exon junctions in c. coordinates (last base of exon before each junction)
        self.junctions_c = []
        if exon_boundaries_c:
            self.junctions_c = sorted(int(x) for x in exon_boundaries_c)
        elif len(self.exons) > 1 and self.coding:
            for s, e in self.exons[:-1]:
                last = e - 1
                k = self.eidx[last]
                self.junctions_c.append(k - self.cds_e0 + 1)
        self.protein_length = len(self.protein.rstrip("*").split("*")[0]) if self.protein else 0

    def _auto_cds(self):
        seq = "".join(self.ref[g] for g in self.exonic)
        prot = translate(seq)
        if seq.startswith("ATG") and len(seq) % 3 == 0 and prot.endswith("*") and "*" not in prot[:-1]:
            s, e = 0, len(seq)
            self.cds_note = "Whole reference used as CDS (ATG..stop, in frame)."
        else:
            orfs = find_orfs(seq, min_aa=10)
            if orfs:
                s, e, _ = orfs[0]
                self.cds_note = f"CDS auto-detected as longest ORF: {s + 1}-{e}."
                self.warnings.append(self.cds_note + " Provide --cds or annotation 'cds' to override.")
            else:
                self.cds_note = "No ORF found."
                self.warnings.append("No ORF found in reference; treating as non-coding.")
                return (0, 0)
        return (self.exonic[s], self.exonic[e - 1] + 1)

    # ------------------------------------------------------------- c. coords
    def _exonic_c(self, g):
        k = self.eidx[g]
        if self.cds_e0 is None:
            return ("n", k + 1)
        if k < self.cds_e0:
            return ("-", self.cds_e0 - k)
        if k > self.cds_e1:
            return ("*", k - self.cds_e1)
        return ("", k - self.cds_e0 + 1)

    @staticmethod
    def _fmt(t):
        kind, n = t
        return f"-{n}" if kind == "-" else (f"*{n}" if kind == "*" else str(n))

    def c_pos(self, g):
        """HGVS c. position string for 0-based reference position g."""
        if g in self.eidx:
            return self._fmt(self._exonic_c(g))
        first, last = self.exons[0][0], self.exons[-1][1] - 1
        if g < first:
            kind, n = self._exonic_c(first)
            d = first - g
            return f"-{n + d}" if kind == "-" else f"-{d}"
        if g > last:
            kind, n = self._exonic_c(last)
            d = g - last
            return f"*{n + d}" if kind == "*" else f"*{d}"
        # intronic
        starts = [s for s, _ in self.exons]
        i = bisect.bisect_right(starts, g) - 1
        up_end = self.exons[i][1] - 1
        dn_start = self.exons[i + 1][0]
        du, dd = g - up_end, dn_start - g
        if du <= dd:
            return f"{self._fmt(self._exonic_c(up_end))}+{du}"
        return f"{self._fmt(self._exonic_c(dn_start))}-{dd}"

    def region(self, g):
        """Return (region, intron_offset) for 0-based position g."""
        if g < 0 or g >= self.L:
            return ("outside", None)
        if g in self.eidx:
            if g in self.cidx:
                return ("cds", 0)
            kind, _ = self._exonic_c(g)
            return ("5utr" if kind == "-" else "3utr" if kind == "*" else "noncoding", 0)
        if g < self.exons[0][0]:
            return ("upstream", None)
        if g >= self.exons[-1][1]:
            return ("downstream", None)
        starts = [s for s, _ in self.exons]
        i = bisect.bisect_right(starts, g) - 1
        du = g - (self.exons[i][1] - 1)
        dd = self.exons[i + 1][0] - g
        return ("intron", du if du <= dd else -dd)

    def near_exon_junction(self, g, window=3):
        """True if an exonic position is within ``window`` bases of an internal exon junction."""
        if g not in self.eidx:
            return False
        if len(self.exons) > 1:
            for s, e in self.exons:
                if not (s <= g < e):
                    continue
                if (s != self.exons[0][0] and g - s < window) or (e != self.exons[-1][1] and e - 1 - g < window):
                    return True
            return False
        if self.junctions_c and g in self.cidx:
            c = self.cidx[g] + 1
            return any(b - window < c <= b + window for b in self.junctions_c)
        return False

    def exon_number(self, g):
        if len(self.exons) == 1 and self.junctions_c and g in self.cidx:
            c = self.cidx[g] + 1
            return bisect.bisect_left(self.junctions_c, c) + 1
        for i, (s, e) in enumerate(self.exons):
            if s <= g < e:
                return i + 1
        if self.junctions_c and g in self.cidx:
            c = self.cidx[g] + 1
            return bisect.bisect_left(self.junctions_c, c) + 1
        return None

    def n_exons(self):
        return len(self.exons) if len(self.exons) > 1 else (len(self.junctions_c) + 1)


# ============================================================ protein effect
def protein_effect(ref_p, alt_p, frameshift=False, hint_pos=None, expected_len_change=0):
    """Compare reference/alternate protein strings (each ending at first '*').

    Returns dict(consequence, hgvs_p, aa_pos, ref_aa, alt_aa, ptc_pos, protein_length).
    """
    n_ref = len(ref_p.rstrip("*"))
    res = dict(consequence=None, hgvs_p="p.?", aa_pos=None, ref_aa="", alt_aa="", ptc_pos=None,
               protein_length=n_ref, aa_change="")
    if ref_p == alt_p:
        p = hint_pos if hint_pos is not None and hint_pos < len(ref_p) else 0
        aa = ref_p[p] if ref_p else "X"
        res.update(consequence="stop_retained_variant" if aa == "*" else "synonymous_variant",
                   hgvs_p=f"p.({a3(aa)}{p + 1}=)", aa_pos=p + 1, ref_aa=aa, alt_aa=aa,
                   aa_change=f"{aa}{p + 1}=")
        return res
    i, L = 0, min(len(ref_p), len(alt_p))
    while i < L and ref_p[i] == alt_p[i]:
        i += 1
    res["aa_pos"] = i + 1
    if i >= len(ref_p):
        res.update(consequence="protein_altering_variant", hgvs_p="p.?")
        return res
    r = ref_p[i]
    res["ref_aa"] = r
    if i == 0 and r == "M":
        res.update(consequence="start_lost", hgvs_p="p.(Met1?)", aa_change="M1?")
        return res
    first_alt = alt_p[i] if i < len(alt_p) else ""
    res["alt_aa"] = first_alt
    if frameshift:
        if first_alt == "*":
            res.update(hgvs_p=f"p.({a3(r)}{i + 1}Ter)", ptc_pos=i + 1, aa_change=f"{r}{i + 1}*")
        elif first_alt == "":
            res.update(hgvs_p=f"p.({a3(r)}{i + 1}fs)")
        else:
            if alt_p.endswith("*"):
                ter = f"Ter{len(alt_p) - i}"
                res["ptc_pos"] = len(alt_p)
            else:
                ter = "Ter?"
            res.update(hgvs_p=f"p.({a3(r)}{i + 1}{a3(first_alt)}fs{ter})",
                       aa_change=f"{r}{i + 1}{first_alt}fs")
        res["consequence"] = "frameshift_variant"
        return res
    if r == "*":
        if alt_p.endswith("*"):
            ext = f"extTer{len(alt_p) - 1 - i}"
        else:
            ext = "extTer?"
        res.update(consequence="stop_lost", hgvs_p=f"p.(Ter{i + 1}{a3(first_alt)}{ext})",
                   aa_change=f"*{i + 1}{first_alt}ext")
        return res
    if first_alt == "*":
        res.update(consequence="stop_gained", hgvs_p=f"p.({a3(r)}{i + 1}Ter)", ptc_pos=i + 1,
                   aa_change=f"{r}{i + 1}*")
        return res
    # premature stop created a few residues downstream by an in-frame delins
    exp_aa = expected_len_change // 3 if expected_len_change % 3 == 0 else 0
    if alt_p.endswith("*") and (len(alt_p) - len(ref_p)) < exp_aa - 1 and len(alt_p) - 1 <= i + 3:
        k = len(alt_p) - 1
        ins = "".join(a3(c) for c in alt_p[i:k])
        res.update(consequence="stop_gained", ptc_pos=k + 1,
                   hgvs_p=f"p.({a3(r)}{i + 1}_{a3(ref_p[k])}{k + 1}delins{ins}Ter)" if k < len(ref_p) else "p.?")
        return res
    j = 0
    while j < len(ref_p) - i and j < len(alt_p) - i and ref_p[-1 - j] == alt_p[-1 - j]:
        j += 1
    D, I = ref_p[i:len(ref_p) - j], alt_p[i:len(alt_p) - j]
    if len(D) == 1 and len(I) == 1:
        res.update(consequence="missense_variant", hgvs_p=f"p.({a3(D)}{i + 1}{a3(I)})", alt_aa=I,
                   aa_change=f"{D}{i + 1}{I}")
    elif not I:
        hg = f"p.({a3(D)}{i + 1}del)" if len(D) == 1 else f"p.({a3(D[0])}{i + 1}_{a3(D[-1])}{i + len(D)}del)"
        res.update(consequence="inframe_deletion", hgvs_p=hg, alt_aa="", aa_change=f"{D}{i + 1}del")
    elif not D:
        k = len(I)
        if i - k >= 0 and ref_p[i - k:i] == I:
            hg = (f"p.({a3(I)}{i}dup)" if k == 1 else f"p.({a3(I[0])}{i - k + 1}_{a3(I[-1])}{i}dup)")
            res["aa_pos"] = i - k + 1
        else:
            hg = f"p.({a3(ref_p[i - 1])}{i}_{a3(ref_p[i])}{i + 1}ins{''.join(a3(c) for c in I)})"
        res.update(consequence="inframe_insertion", hgvs_p=hg, alt_aa=I, ref_aa="", aa_change=f"{i}ins{I}")
    else:
        rng = f"{a3(D[0])}{i + 1}" if len(D) == 1 else f"{a3(D[0])}{i + 1}_{a3(D[-1])}{i + len(D)}"
        res.update(consequence="protein_altering_variant",
                   hgvs_p=f"p.({rng}delins{''.join(a3(c) for c in I)})", alt_aa=I, ref_aa=D,
                   aa_change=f"{D}{i + 1}delins{I}")
    return res


# ============================================================ HGVS helpers
def hgvs_c(v, tm: TranscriptModel, prefix="c."):
    cp = tm.c_pos
    s, e = v.ref_start, v.ref_end
    if v.kind == "SNV":
        return f"{prefix}{cp(s)}{v.ref}>{v.alt}"
    if v.kind == "MNV" or v.kind == "delins":
        rng = cp(s) if e - s == 1 else f"{cp(s)}_{cp(e - 1)}"
        return f"{prefix}{rng}delins{v.alt}"
    if v.kind == "deletion":
        rng = cp(s) if e - s == 1 else f"{cp(s)}_{cp(e - 1)}"
        return f"{prefix}{rng}del"
    if v.kind == "duplication":
        k = len(v.alt)
        rng = cp(s - 1) if k == 1 else f"{cp(s - k)}_{cp(s - 1)}"
        return f"{prefix}{rng}dup"
    if v.kind == "insertion":
        return f"{prefix}{cp(s - 1)}_{cp(s)}ins{v.alt}"
    return f"{prefix}?"


def hgvs_raw(v, prefix="g."):
    s, e = v.ref_start + 1, v.ref_end
    if v.kind == "SNV":
        return f"{prefix}{s}{v.ref}>{v.alt}"
    if v.kind in ("MNV", "delins"):
        return f"{prefix}{s}{'' if e == s else '_' + str(e)}delins{v.alt}"
    if v.kind == "deletion":
        return f"{prefix}{s}{'' if e == s else '_' + str(e)}del"
    if v.kind == "duplication":
        k = len(v.alt)
        return f"{prefix}{s - k}{'' if k == 1 else '_' + str(s - 1)}dup"
    if v.kind == "insertion":
        return f"{prefix}{s - 1}_{s}ins{v.alt}"
    return prefix + "?"


def hgvs_p_protein_mode(v):
    """HGVS p. for variants called directly on protein sequences."""
    return None  # computed through protein_effect


# ============================================================ annotation
def annotate_dna(v, tm: TranscriptModel):
    """Return dict with consequence info for a DNA RawVariant."""
    out = {"hgvs_c": hgvs_c(v, tm), "hgvs_ref": hgvs_raw(v), "region": None, "intron_offset": None,
           "consequence": "non_coding_variant", "hgvs_p": "", "aa_pos": None, "ref_aa": "", "alt_aa": "",
           "ptc_pos": None, "protein_length": tm.protein_length, "codon_change": "", "exon": None,
           "splice_region": False, "aa_change": ""}
    s, e = v.ref_start, v.ref_end
    if v.kind in ("insertion", "duplication"):
        positions = [p for p in (s - 1, s) if 0 <= p < tm.L]
    else:
        positions = list(range(s, e))
    regions = [tm.region(g) for g in positions]
    out["exon"] = tm.exon_number(positions[0]) if positions else None
    # splice assessment
    canon, splice_reg = None, False
    for (rg, off) in regions:
        if rg == "intron" and off is not None:
            if abs(off) <= 2:
                canon = "splice_donor_variant" if off > 0 else "splice_acceptor_variant"
            elif abs(off) <= 8:
                splice_reg = True
    if any(tm.near_exon_junction(g) for g in positions):
        splice_reg = True
    out["splice_region"] = splice_reg or bool(canon)
    kinds = {rg for rg, _ in regions}
    first_reg = regions[0][0] if regions else "outside"
    out["region"] = "/".join(sorted(kinds))
    intr = [off for rg, off in regions if rg == "intron"]
    out["intron_offset"] = intr[0] if intr else None

    all_coding = positions and all(rg == "cds" for rg, _ in regions)
    if v.kind in ("insertion", "duplication") and positions and positions[-1] in tm.cidx and len(positions) == 1:
        all_coding = True
    if all_coding and v.kind not in ("insertion", "duplication"):
        ci = tm.cidx[s]
        if tm.cidx[e - 1] - ci != e - 1 - s:
            all_coding = False  # spans an intron
    if all_coding:
        if v.kind in ("insertion", "duplication"):
            ci = cj = tm.cidx[s] if s in tm.cidx else tm.cidx[s - 1] + 1
        else:
            ci, cj = tm.cidx[s], tm.cidx[e - 1] + 1
        alt_cds = tm.cds_seq[:ci] + v.alt + tm.cds_seq[cj:]
        ref_p = translate(tm.cds_seq + tm.utr3_seq, to_stop=True)
        alt_p = translate(alt_cds + tm.utr3_seq, to_stop=True)
        fs = (len(v.alt) - (cj - ci)) % 3 != 0
        pe = protein_effect(ref_p, alt_p, frameshift=fs, hint_pos=ci // 3,
                            expected_len_change=len(v.alt) - (cj - ci))
        out.update({k: pe[k] for k in ("consequence", "hgvs_p", "aa_pos", "ref_aa", "alt_aa", "ptc_pos", "aa_change")})
        out["protein_length"] = pe["protein_length"]
        if v.kind == "SNV":
            c0 = (ci // 3) * 3
            rc = tm.cds_seq[c0:c0 + 3]
            ac = alt_cds[c0:c0 + 3]
            out["codon_change"] = f"{rc}>{ac}"
        if canon:
            out["consequence"] = canon
        elif splice_reg and out["consequence"] in ("synonymous_variant", "missense_variant"):
            out["secondary_consequence"] = "splice_region_variant"
        return out
    if canon:
        out["consequence"] = canon
        return out
    coding_hit = [g for g in positions if g in tm.cidx]
    if v.kind in ("insertion", "duplication") and len(coding_hit) == 1 and "intron" not in kinds:
        # insertion right at the CDS edge: lies in the UTR
        out["consequence"] = "5_prime_UTR_variant" if coding_hit[0] == positions[-1] else "3_prime_UTR_variant"
        return out
    if coding_hit:
        # variant partially overlaps CDS (spans exon/intron or CDS boundary)
        if any(rg == "intron" for rg, _ in regions):
            out["consequence"] = "exon_loss_variant"
        elif tm.cidx[coding_hit[0]] < 3:
            out["consequence"] = "start_lost"
            out["hgvs_p"] = "p.(Met1?)"
        else:
            out["consequence"] = "stop_lost"
            out["hgvs_p"] = "p.?"
        return out
    mapping = {"5utr": "5_prime_UTR_variant", "3utr": "3_prime_UTR_variant", "intron": "intron_variant",
               "upstream": "upstream_gene_variant", "downstream": "downstream_gene_variant"}
    out["consequence"] = "splice_region_variant" if splice_reg else mapping.get(first_reg, "non_coding_variant")
    return out


def annotate_protein(v, ref_prot):
    """Consequence for a variant called on a protein-vs-protein alignment."""
    ref_p = ref_prot.split("*")[0] + "*"
    alt_full = ref_prot[:v.ref_start] + v.alt + ref_prot[v.ref_end:]
    open_end = v.kind in ("frameshift", "extension") and "*" not in v.alt
    alt_p = alt_full if open_end else alt_full.split("*")[0] + "*"
    pe = protein_effect(ref_p, alt_p, frameshift=(v.kind == "frameshift"), hint_pos=v.ref_start,
                        expected_len_change=3 * (len(v.alt) - len(v.ref)))
    out = {"hgvs_c": "", "hgvs_ref": "", "region": "protein", "intron_offset": None,
           "codon_change": "", "exon": None, "splice_region": False}
    out.update(pe)
    if v.kind in ("frameshift", "extension") and pe.get("aa_pos"):
        out["alt_tail"] = alt_p[pe["aa_pos"] - 1:]
    return out
