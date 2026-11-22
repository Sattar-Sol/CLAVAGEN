CELF2 test family FAM-CELF2-01 - sequencing files for CLAVAGEN
=====================================================================

Reference: NM_001025077.3 (CELF2 transcript variant 3, protein NP_001020248.1, 508 aa),
taken from NCBI on 2026-11-21 and checked (identical SHA-256; translation = NCBI protein).
The patient sequences are SIMULATED from this reference (cDNA amplicons and synthetic traces).

Who is who
----------
Member      Relationship  Status      Genotype (NM_001025077.3)
P1_Father   father        unaffected  c.1478G>A p.(Arg493His)                het
P2_Mother   mother        unaffected  c.688C>T  p.(Gln230Ter)                het
P3_Child1   sibling       unaffected  no variant
P4_Child2   sibling       unaffected  c.1478G>A p.(Arg493His)                het (from father)
P5_Child3   sibling       unaffected  c.688C>T  p.(Gln230Ter)                het (from mother)
P6_Child4   PROBAND       AFFECTED    c.688C>T + c.1478G>A, in trans         compound het
P7_Child5   sibling       unaffected  c.1141del p.(Asp381ThrfsTer82)         het, de novo

Same variants on the MANE Select transcript NM_001326342.2 (NP_001313271.1):
  c.1478G>A p.(Arg493His)  = NM_001326342.2 p.(Arg506His) (c.1517G>A)
  c.688C>T  p.(Gln230Ter)  = NM_001326342.2 p.(Gln237Ter) (c.709C>T)
The software names variants on the reference you load (here NM_001025077.3).

Files (formats are mixed on purpose)
------------------------------------
0_Reference/  CELF2_NM_001025077.3.fasta          DNA reference (nt 1-2203: 5'UTR, CDS 77-1603, 3'UTR start)
              CELF2_NP_001020248.1_protein.fasta  protein reference
              CELF2_annotation_AR_optional.json   optional: inheritance AR + CDS + exon junctions
1_DNA/<patient>/
  P1  .ab1 amplicon A forward + amplicon B REVERSE
  P2  .fasta  full-CDS consensus
  P3  .seq    plain text, no header
  P4  .ab1 amplicons A + B forward
  P5  .fa     multi-FASTA, two amplicon consensuses
  P6  .ab1 A + B forward and a .txt consensus
  P7  .ab1 A forward, B forward AND B reverse (heterozygous deletion: double peaks)
  Amplicon A = c.501-1000 (covers c.688), amplicon B = c.1001-1560 (covers c.1141 and c.1478).
  Consensus files use IUPAC codes for heterozygous positions (Y = C/T at c.688, R = A/G at c.1478).
2_Protein/<patient>_CELF2_protein.fasta  one protein per allele; truncated proteins end with '*'.

Patient lists (CSV) - family, relationship, affected status already filled in
------------------------------------------------------------------------------
CELF2_DNA_patients_IMPORT_in_window.csv      DNA, sequences written inside the CSV (window import).
    Each person's two amplicons as IUPAC consensus. A heterozygous deletion cannot be written as
    text, so P7 shows c.1141del as "apparent homozygous": after import, Edit P7, clear amplicon B
    from the sequence box and attach 1_DNA/P7_Child5/P7_Child5_CELF2_ampB_F.ab1 (and _R.ab1).
CELF2_DNA_patients_command_line.csv          DNA, paths to the real files incl. .ab1 traces
    (command line; paths are relative to the CSV). Run from the CLAVAGEN folder:
    python3 -m clavagen cohort --patients examples/CELF2_family/CELF2_DNA_patients_command_line.csv \
        -r examples/CELF2_family/0_Reference/CELF2_NM_001025077.3.fasta \
        --ref-setting CELF2:inheritance=AR -o results/
CELF2_PROTEIN_patients_IMPORT_in_window.csv  protein, one protein per allele inside the CSV
    (works in the window and on the command line).
The window cannot open files named inside a CSV (it only receives the CSV's text) - that is why
the window lists carry the sequences themselves.

Window, with a CSV: 1. Reference > add 0_Reference/CELF2_NM_001025077.3.fasta (Inheritance AR)
  4. Patients > Import patients > choose the CSV > Run. (Protein: same, protein CSV; either reference.)

How to run in the window (adding patients by hand)
--------------------------------------------------
1. Reference: Add FASTA file(s) > 0_Reference/CELF2_NM_001025077.3.fasta; set Inheritance to AR
   (or load the optional annotation file in section 3).
2. Keep "Search online databases automatically" on.
3. + Add patient for each member: ID, Family ID FAM-CELF2-01, relationship and affected status
   as in the table (P6 = proband). Choose files: FASTA/.seq/.txt/.fa fill the sequence box,
   .ab1 files are attached as traces.
4. Run. Protein test: same steps with 2_Protein files (against either reference; the DNA
   reference lets the DNA change be inferred offline).

Expected (tested with CLAVAGEN 1.0.0)
--------------------------
DNA: all carriers and the proband called heterozygous from every file type; P6 in trans
(both parents tested); P7 c.1141del heterozygous from the forward AND the reverse trace and
"de novo in sibling"; P3 no variants. Protein: p.(Gln230Ter) from the '*', the frameshift
p.(Asp381ThrfsTer82) with c.1141del inferred, zygosity n.d.
Classes and couple risk need online annotation (ClinVar/gnomAD); offline they stay VUS.

Note: CELF2 disease (DEE97) is autosomal dominant; this recessive family is a constructed test.
p.(Gln230Ter) and c.1141del are invented test variants.
