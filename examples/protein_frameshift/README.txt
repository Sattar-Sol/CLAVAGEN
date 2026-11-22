Protein-level truncating variants (V11)

Reference: HBB_reference_NM_000518.5_CDS.fa (DNA, so the DNA change can be inferred offline)
Query:     patients_protein_frameshift.fa (two protein sequences, each ending with *)

Expected:
  Patient_FS        p.(Phe42LeufsTer19)  frameshift, inferred c.126_129del, PVS1
  Patient_Nonsense  p.(Gln40Ter)         premature stop, inferred c.118C>T, PVS1

Command line:
  python3 -m clavagen run -q patients_protein_frameshift.fa -r HBB_reference_NM_000518.5_CDS.fa -g HBB -o results/
In the window: load the reference in section 1 and paste the query file into "Add patient".
End a complete protein with * - without it a short query is treated as a partial sequence.
