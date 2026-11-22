"""Sequence input: FASTA files, plain sequence files and literal sequences."""

from __future__ import annotations

import hashlib
import os
import re
from dataclasses import dataclass, field

IUPAC_DNA = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
# Two-allele codes -> heterozygous genotype in Sanger / consensus sequences
HET_CODES = {k: v for k, v in IUPAC_DNA.items() if len(v) == 2}
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBZXUO*")


@dataclass
class SeqRecord:
    id: str
    seq: str
    description: str = ""
    mol_type: str = "dna"  # dna | protein
    source: str = ""
    sha256: str = ""
    warnings: list = field(default_factory=list)

    def __len__(self):
        return len(self.seq)


def parse_fasta_text(text: str, default_id: str = "seq"):
    records, cur, desc, buf = [], None, "", []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(";") or line.startswith("#"):
            continue
        if line.startswith(">"):
            if cur is not None:
                records.append((cur, desc, "".join(buf)))
            header = line[1:].strip()
            parts = header.split(None, 1)
            cur = parts[0] if parts else f"{default_id}{len(records) + 1}"
            desc = parts[1] if len(parts) > 1 else ""
            buf = []
        else:
            if cur is None:
                cur = default_id
            buf.append(re.sub(r"[\s\d/]", "", line))
    if cur is not None:
        records.append((cur, desc, "".join(buf)))
    return records


def detect_mol_type(seq: str) -> str:
    s = [c for c in seq.upper() if c.isalpha()]
    if not s:
        return "dna"
    acgtn = sum(1 for c in s if c in "ACGTUN")
    all_iupac = all(c in IUPAC_DNA or c == "U" for c in s)
    return "dna" if (all_iupac and acgtn / len(s) >= 0.85) else "protein"


def clean_sequence(seq: str, mol_type: str):
    seq = re.sub(r"[\s\d\-\.]", "", seq).upper()
    warnings = []
    if mol_type == "dna":
        if "U" in seq:
            seq = seq.replace("U", "T")
            warnings.append("RNA input: U converted to T")
        bad = sorted({c for c in seq if c not in IUPAC_DNA})
        if bad:
            warnings.append(f"Invalid DNA characters replaced with N: {''.join(bad)}")
            seq = "".join(c if c in IUPAC_DNA else "N" for c in seq)
        n_amb = sum(1 for c in seq if c not in "ACGT")
        if n_amb:
            warnings.append(f"{n_amb} ambiguous/IUPAC base(s) present")
    else:
        bad = sorted({c for c in seq if c not in PROTEIN_ALPHABET})
        if bad:
            warnings.append(f"Invalid protein characters replaced with X: {''.join(bad)}")
            seq = "".join(c if c in PROTEIN_ALPHABET else "X" for c in seq)
        # keep internal content up to the first stop as the reference protein
    return seq, warnings


def load_sequences(specs, role: str = "query", force_type: str | None = None):
    """Load sequences from a list of specs.

    Each spec may be: a FASTA/plain-text file path, a literal ``NAME=SEQUENCE``
    or a bare literal sequence. Returns a list of SeqRecord.
    """
    out = []
    for spec in specs or []:
        spec = spec.strip()
        if not spec:
            continue
        if os.path.isfile(spec):
            with open(spec, "rb") as fh:
                head = fh.read(4)
            if head == b"ABIF":
                from .trace import load_trace
                with open(spec, "rb") as fh:
                    tr = load_trace(fh.read(), os.path.splitext(os.path.basename(spec))[0])
                seq = tr.consensus()
                rec = SeqRecord(id=os.path.splitext(os.path.basename(spec))[0], seq=seq,
                                description=f"Sanger trace {os.path.basename(spec)} ({tr.name})", mol_type="dna",
                                source=os.path.abspath(spec), sha256=hashlib.sha256(open(spec, "rb").read()).hexdigest(),
                                warnings=list(tr.notes) + [f"Trace trimmed to peaks {tr.trim[0] + 1}-{tr.trim[1]} of "
                                                           f"{len(tr.peaks)}; {sum(tr.mixed[tr.trim[0]:tr.trim[1]])} double "
                                                           "peak(s) in the kept region."])
                rec.trace = tr
                out.append(rec)
                continue
            with open(spec, "r", encoding="utf-8", errors="replace") as fh:
                text = fh.read()
            base = os.path.splitext(os.path.basename(spec))[0]
            recs = parse_fasta_text(text, default_id=base)
            src = os.path.abspath(spec)
        else:
            if spec.startswith(">"):
                recs = parse_fasta_text(spec.replace("\\n", "\n"), default_id=f"{role}{len(out) + 1}")
            elif "=" in spec and re.match(r"^[\w.\-|:]+=", spec):
                name, s = spec.split("=", 1)
                recs = [(name, "", s)]
            else:
                recs = [(f"{role}{len(out) + 1}", "", spec)]
            src = "literal"
        for rid, desc, s in recs:
            mol = force_type if force_type in ("dna", "protein") else detect_mol_type(s)
            seq, warns = clean_sequence(s, mol)
            if not seq:
                continue
            out.append(SeqRecord(
                id=rid, seq=seq, description=desc, mol_type=mol, source=src,
                sha256=hashlib.sha256(seq.encode()).hexdigest(), warnings=warns,
            ))
    return out


def write_fasta(records, path, width=60):
    with open(path, "w") as fh:
        for r in records:
            fh.write(f">{r.id} {r.description}".rstrip() + "\n")
            for i in range(0, len(r.seq), width):
                fh.write(r.seq[i:i + width] + "\n")
