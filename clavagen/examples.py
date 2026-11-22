"""Worked demonstration dataset based on the human HBB coding sequence (NM_000518.5 CDS).

Population frequencies, in-silico scores and phenotype inputs in the demo
annotation are ILLUSTRATIVE ONLY and must not be used clinically.
"""

from __future__ import annotations

import json
import os

HBB_CDS = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGG"
    "CTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAG"
    "GCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCAC"
    "TGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACC"
    "CCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
)
# Human SOD1 coding sequence (translation identical to NP_000445.1)
SOD1_CDS = ("ATGGCGACGAAGGCCGTGTGCGTGCTGAAGGGCGACGGCCCAGTGCAGGGCATCATCAATTTCGAGCAGAAGGAAAGTAATGGACCAGTGAAGGTGTGGGGAAGC"
       "ATTAAAGGACTGACTGAAGGCCTGCATGGATTCCATGTTCATGAGTTTGGAGATAATACAGCAGGCTGTACCAGTGCAGGTCCTCACTTTAATCCTCTATCCAG"
       "AAAACACGGTGGGCCAAAGGATGAAGAGAGGCATGTTGGAGACTTGGGCAATGTGACTGCTGACAAAGATGGTGTGGCCGATGTGTCTATTGAAGATTCTGTG"
       "ATCTCACTCTCAGGAGACCATTGCATCATTGGCCGCACACTGGTGGTCCATGAAAAAGCAGATGACTTGGGCAAAGGTGGAAATGAAGAAAGTACAAAGACAG"
       "GAAACGCTGGAAGTCGTTTGGCTTGTGGTGTAATTGGGATCGCCCAATAA")
HBB_PROTEIN = ("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATL"
               "SELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH")


def _mut(seq, pos1, ref, alt):
    i = pos1 - 1
    assert seq[i:i + len(ref)] == ref, (pos1, ref, seq[i:i + len(ref)])
    return seq[:i] + alt + seq[i + len(ref):]


def _wrap(s, w=60):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


def write_examples(outdir):
    os.makedirs(outdir, exist_ok=True)
    # Patient 1: heterozygous c.20A>T (IUPAC W = A/T) - HbS carrier
    p1 = _mut(HBB_CDS, 20, "A", "W")
    # Patient 2: homozygous c.118C>T p.Gln40Ter (beta0) + synonymous c.9T>C
    p2 = _mut(_mut(HBB_CDS, 118, "C", "T"), 9, "T", "C")
    # Patient 3: partial amplicon (c.51-c.360) carrying c.126_129del (legacy codons 41/42 -TTCT)
    p3 = _mut(HBB_CDS, 124, "TTCT", "")[50:356]
    # Patient 4: in-frame deletion of Lys at codon 18 (c.52_54del) + low-quality N base
    p4 = _mut(HBB_CDS, 52, "AAG", "")
    p4 = p4[:296] + "N" + p4[297:]
    ref_p = os.path.join(outdir, "HBB_reference.fa")
    with open(ref_p, "w") as fh:
        fh.write(f">HBB_NM_000518.5 Homo sapiens hemoglobin subunit beta CDS\n{_wrap(HBB_CDS)}\n")
    q_p = os.path.join(outdir, "patients_dna.fa")
    with open(q_p, "w") as fh:
        fh.write(f">Patient1_HBB Sanger consensus (heterozygous calls as IUPAC)\n{_wrap(p1)}\n")
        fh.write(f">Patient2_HBB\n{_wrap(p2)}\n")
        fh.write(f">Patient3_HBB_amplicon partial amplicon\n{_wrap(p3)}\n")
        fh.write(f">Patient4_HBB\n{_wrap(p4)}\n")
    prot_p = os.path.join(outdir, "patient_protein.fa")
    with open(prot_p, "w") as fh:
        fh.write(f">Patient5_HBB_protein mass-spec/derived protein\n{_wrap(HBB_PROTEIN[:6] + 'V' + HBB_PROTEIN[7:])}\n")
    ann = {
        "_comment": "DEMONSTRATION annotation. Population/in-silico values are illustrative, not real data.",
        "gene": "HBB", "transcript": "NM_000518.5", "protein_id": "NP_000509.1",
        "chrom": "HBB_NM_000518.5", "cds": [1, 444], "exon_boundaries_c": [92, 315],
        "inheritance": "AR", "lof_mechanism": True,
        "domains": [{"name": "Globin", "start": 3, "end": 147, "critical": False},
                    {"name": "Heme pocket (illustrative)", "start": 63, "end": 92, "critical": True}],
        "known_variants": [
            {"hgvs_c": "c.20A>T", "hgvs_p": "p.Glu7Val", "classification": "Pathogenic", "source": "ClinVar (HbS)"},
            {"hgvs_c": "c.118C>T", "hgvs_p": "p.Gln40Ter", "classification": "Pathogenic", "source": "ClinVar (beta0)"},
            {"hgvs_c": "c.19G>A", "hgvs_p": "p.Glu7Lys", "classification": "Pathogenic", "source": "ClinVar (HbC)"},
        ],
        "population_dataset": "gnomAD (illustrative values)",
        "population": {"c.20A>T": {"af": 0.004, "source": "illustrative"},
                       "c.9T>C": {"af": 0.02, "source": "illustrative"}},
        "ba1_exceptions": ["c.20A>T"],
        "insilico": {},
    }
    ann_p = os.path.join(outdir, "HBB_annotation.json")
    with open(ann_p, "w") as fh:
        json.dump(ann, fh, indent=2)
    return {"reference_dna": ref_p, "dna_queries": q_p, "protein_query": prot_p, "annotation": ann_p}
