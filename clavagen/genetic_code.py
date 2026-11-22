"""Genetic code, amino-acid properties and substitution matrices."""

from __future__ import annotations

import math

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AAS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter",
    "X": "Xaa", "U": "Sec", "B": "Asx", "Z": "Glx", "O": "Pyl",
}
AA1 = {v: k for k, v in AA3.items()}
AA_NAME = {
    "A": "alanine", "R": "arginine", "N": "asparagine", "D": "aspartate", "C": "cysteine",
    "Q": "glutamine", "E": "glutamate", "G": "glycine", "H": "histidine", "I": "isoleucine",
    "L": "leucine", "K": "lysine", "M": "methionine", "F": "phenylalanine", "P": "proline",
    "S": "serine", "T": "threonine", "W": "tryptophan", "Y": "tyrosine", "V": "valine",
    "*": "stop",
}

# Physico-chemical class used to describe the nature of a substitution.
AA_CLASS = {
    "A": "nonpolar", "V": "nonpolar", "L": "nonpolar", "I": "nonpolar", "M": "nonpolar",
    "F": "aromatic", "W": "aromatic", "Y": "aromatic",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar",
    "K": "positive", "R": "positive", "H": "positive",
    "D": "negative", "E": "negative",
    "G": "special", "P": "special", "C": "special",
    "*": "stop",
}
# Kyte-Doolittle hydropathy
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

_BLOSUM62_TXT = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""


def _parse_blosum(txt: str) -> dict:
    lines = [ln for ln in txt.strip("\n").splitlines() if ln.strip()]
    cols = lines[0].split()
    mat = {}
    for ln in lines[1:]:
        parts = ln.split()
        r = parts[0]
        for c, v in zip(cols, parts[1:]):
            mat[(r, c)] = int(v)
    return mat


BLOSUM62 = _parse_blosum(_BLOSUM62_TXT)


def blosum62(a: str, b: str) -> int:
    a = a if (a, a) in BLOSUM62 else "X"
    b = b if (b, b) in BLOSUM62 else "X"
    return BLOSUM62[(a, b)]


# Grantham (1974) physico-chemical distance, computed from the original
# composition (c), polarity (p) and volume (v) values with the published
# weights and scaling factor. Values match the published table within ~1-2 units.
_GRANTHAM_PROPS = {
    "S": (1.42, 9.2, 32.0), "R": (0.65, 10.5, 124.0), "L": (0.0, 4.9, 111.0),
    "P": (0.39, 8.0, 32.5), "T": (0.71, 8.6, 61.0), "A": (0.0, 8.1, 31.0),
    "V": (0.0, 5.9, 84.0), "G": (0.74, 9.0, 3.0), "I": (0.0, 5.2, 111.0),
    "F": (0.0, 5.2, 132.0), "Y": (0.20, 6.2, 136.0), "C": (2.75, 5.5, 55.0),
    "H": (0.58, 10.4, 96.0), "Q": (0.89, 10.5, 85.0), "N": (1.33, 11.6, 56.0),
    "K": (0.33, 11.3, 119.0), "D": (1.38, 13.0, 54.0), "E": (0.92, 12.3, 83.0),
    "M": (0.0, 5.7, 105.0), "W": (0.13, 5.4, 170.0),
}


def grantham(a: str, b: str):
    """Return Grantham distance (0-215) or None for non-standard residues."""
    if a not in _GRANTHAM_PROPS or b not in _GRANTHAM_PROPS:
        return None
    c1, p1, v1 = _GRANTHAM_PROPS[a]
    c2, p2, v2 = _GRANTHAM_PROPS[b]
    d = math.sqrt(1.833 * (c1 - c2) ** 2 + 0.1018 * (p1 - p2) ** 2 + 0.000399 * (v1 - v2) ** 2)
    return int(round(50.723 * d))


def grantham_category(d):
    if d is None:
        return "n/a"
    if d <= 50:
        return "conservative"
    if d <= 100:
        return "moderately conservative"
    if d <= 150:
        return "moderately radical"
    return "radical"


COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn-", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn-")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq: str, to_stop: bool = False) -> str:
    """Translate a DNA sequence (frame 0). Incomplete trailing codons are ignored.

    If ``to_stop`` the translation is truncated *after* the first stop codon
    (the '*' is kept so callers can see that a stop was reached).
    """
    seq = seq.upper().replace("U", "T")
    out = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON_TABLE.get(seq[i:i + 3], "X")
        out.append(aa)
        if to_stop and aa == "*":
            break
    return "".join(out)


def find_orfs(seq: str, min_aa: int = 30):
    """Forward-strand ORFs (ATG..stop). Returns list of (start0, end0_exclusive, protein)."""
    seq = seq.upper()
    orfs = []
    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            if seq[i:i + 3] == "ATG":
                j = i
                while j < len(seq) - 2:
                    if CODON_TABLE.get(seq[j:j + 3]) == "*":
                        break
                    j += 3
                if j < len(seq) - 2:
                    prot = translate(seq[i:j + 3])
                    if len(prot) - 1 >= min_aa:
                        orfs.append((i, j + 3, prot))
                    i = j + 3
                    continue
                break
            i += 3
    orfs.sort(key=lambda o: o[1] - o[0], reverse=True)
    return orfs


def gc_content(seq: str) -> float:
    s = [c for c in seq.upper() if c in "ACGT"]
    if not s:
        return 0.0
    return 100.0 * sum(1 for c in s if c in "GC") / len(s)


def three(aa_seq: str) -> str:
    return "".join(AA3.get(a, "Xaa") for a in aa_seq)
