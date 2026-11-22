# Changelog

## CLAVAGEN 1.0.0 (2026-11-21) - first public release
The project is renamed **CLAVAGEN** (CLassification, Alignment, VAriant calling, GENotype & GENealogy).
It was formerly **ACMG-VarCall**; this release is functionally ACMG-VarCall V12.

- **New name everywhere:**
  - package `clavagen` and command `python3 -m clavagen ...`;
  - app `CLAVAGEN.app`, window title and report headers;
  - default results folder `Documents/CLAVAGEN Results`.
- **Old name still works:** `python3 -m acmg_varcall ...` and `import acmg_varcall` forward to CLAVAGEN.
- **Settings carried over:** on first start, `~/.acmg_varcall` (settings, reference library, database cache) is copied to `~/.clavagen`. Your saved results folder and history are kept.
- **GitHub files:** LICENSE, `.gitignore`, `CITATION.cff` and a GitHub README with screenshots.
- **Demo results** in `examples/` regenerated with this release.
- **CELF2 family** added as a worked example (`examples/CELF2_family/`: mixed file types, compound heterozygous proband, de novo frameshift).
- **Fixes:**
  - **Gene symbol offline:** it is now read from the reference header (`gene=CELF2`, `(CELF2)`, `GN=CELF2`) on the command line too, not only in the window. Without it the genotype step was skipped.
  - **Genotype without a gene symbol:** variants are now grouped per reference sequence instead.
  - **`--ref-setting` by gene symbol:** it accepts a gene symbol as well as a reference ID (e.g. `CELF2:inheritance=AR`).
  - **Indels in a consensus/text sequence:** these are reported with *unresolved* zygosity instead of homozygous, because such a sequence cannot show a heterozygous indel.
  - **Cohort overview:** variants found in several files of one patient are listed once.
  - **Pedigree:** the proband arrow no longer overlaps the name.

The versions below were released as `ACMG-VarCall_V<n>`.

## V12 (2026-11-21)
- **Genotype per patient, not per file.** All of a patient's variants in a gene are combined across their files (amplicons, forward/reverse traces, consensus). Before, a compound heterozygote whose two variants were in different files was never recognised, because each file was interpreted alone.
- **Compound heterozygosity reported whatever the classification.** Two heterozygous variants are described as:
  - *compound heterozygous, in trans* (phase from the family analysis);
  - *in cis*;
  - *phase unknown (possible compound heterozygote)*, with a recommendation to establish phase.
  Before, this was only mentioned when both variants were P/LP. The class of each allele is stated next to it (both P/LP, one P/LP + one VUS, or both VUS), with the matching clinical wording.
- **Variant (mutation) types in the genotype:**
  - each allele is labelled null (nonsense, frameshift, canonical splice, start-loss), missense, in-frame, stop-loss, splice region or synonymous;
  - the genotype states the combination (e.g. *null (nonsense) / missense*);
  - comments cover null/null, null/missense and missense/missense genotypes;
  - a note explains when a null allele could not get PVS1 because loss of function is not established for the gene.
- **New "Genotype" box** in the report (and line in the Markdown summary). The cohort overview status now leads with the genotype, e.g. "Compound heterozygous, in trans (confirmed by testing both parents) - both variants VUS [null (nonsense) / missense]".
- Single heterozygous VUS are shown as "Heterozygous VUS [type]" instead of "No reportable finding".
- Protein input: two variants on separate input sequences are flagged as a possible compound heterozygote (confirm at DNA level).
- Version 12.0.0.

## V11 (2026-11-21)
- **Fixed: frameshifts in protein queries were missed.** A frameshifted protein matches the reference only up to the break, so the local alignment stopped there and the report said "No variants detected" (a false negative). Protein queries are now checked for what the alignment left out:
  - **Frameshift:** the query continues with new residues (and usually stops early), e.g. `p.(Phe42LeufsTer19)`. Residues that aligned by chance just after the break are absorbed into the frameshift rather than called as separate missense variants. Without a stop in the supplied sequence the name ends in `fsTer?`.
  - **Premature stop:** a query ending in `*` before the reference end is called as nonsense, e.g. `p.(Gln40Ter)`.
  - **Stop-loss:** residues beyond the reference stop are called as an extension, e.g. `p.(Ter148GlyextTer?)`.
  - **Partial sequence:** a query that simply ends early gets a QC warning naming the residues not analysed (end a complete protein with `*`).
- **DNA change inferred for protein frameshifts:** every 1, 2, 4 or 5 nt deletion and 1 or 2 nt insertion near the break is tested against the transcript (DNA reference or RefSeq). A unique explanation is reported (e.g. `c.126_129del`, 3'-normalised) and used for database lookups; when codon redundancy allows several (e.g. `c.50dup` or `c.50_51insA`), all are listed.
- Nucleotide inference for protein queries now also works offline when the reference is a DNA sequence.
- Exon junctions from the annotation or RefSeq are used for PVS1/NMD prediction with protein references.
- **Fixed: de novo variants in siblings were not reported.** De novo checks ran only on the proband's variants. A variant in a sibling that is absent in both parents is now reported as de novo; PM6 (or PS2 with confirmed parentage) is applied when that sibling is affected.
- **Fixed: heterozygous indels in reverse Sanger reads.** When a reverse read is turned round, the out-of-register double peaks lie before the indel rather than after it. They were called as dozens of false heterozygous SNVs; the alleles are now separated correctly and one heterozygous indel is called (forward and reverse reads agree).
- Version 11.0.0.

## V10 (2026-11-20)
- **Fixed:** importing a second patient list (e.g. the HBB family after the SOD1 family) added no patients. Patients whose IDs were already in the list (`P1_Father`, `P2_Mother`...) were skipped silently, so it looked like a 7-patient limit. There is no patient limit.
- **New dialog for duplicate IDs on import:**
  - **Rename new patients** (the Family ID is added in front, e.g. `FAM-HBB-01_P1_Father`; recommended);
  - **Replace existing**;
  - **Skip duplicates**;
  - **Cancel**.
  - The message after import lists exactly what was renamed, replaced or skipped.
- **Clashing Family IDs** from a different list are detected and renamed (e.g. `FAM-HBB-01_2`) unless you confirm they are the same family, so separate families are not merged.
- **Several patient lists** can be selected in one import.
- Version 10.0.0.

## V9 (2026-11-20)
- **Several references in one run:** section 1 is now a reference list.
  - Add references from multiple FASTA files, RefSeq accessions or pasted text.
  - The gene is detected per reference.
  - Settings per reference: gene, inheritance (AD, AR, AD/AR, XLR, XLD) and CDS range. No annotation file is needed.
  - Preview works for each reference or for all of them.
  - A duplicate reference ID is reported.
- **Reference library:** save references and their settings for future sessions, and add them back with one click.
- **Patient → gene:** optional "Reference / gene" choice in the pop-up and a `reference` CSV column; auto-matching stays the default.
- **Cohort overview** shows the gene for every variant.
- **Command line:** new `--ref-setting REF_ID:key=value` option for `run` and `cohort` (several `-r` allowed).
- **Per-reference settings** match with or without the version number (NM_000454 = NM_000454.5).
- **"Load multi-gene example":** HBB + SOD1, two families, 12 people.
- **Tests:** 33 automated tests.
- Version 9.0.0.

## V8 (2026-11-20)
- **Family and couple analysis:** Family ID, relationship and affected status per patient (pop-up, CSV).
  - De novo evidence (PS2/PM6).
  - Phase in the proband (PM3 in trans / BP2 in cis), with compound-heterozygote detection.
  - Segregation (PP1 / BS4).
  - Reduced-penetrance notes.
  - Couple reproductive risk (proband + partner, mother + father).
  - Members are re-classified with the family evidence.
  - Outputs: pedigree drawing, segregation table, family section in reports, "Load family example".
- **Sanger trace (.ab1) input:**
  - ABIF reader, verified against Biopython on real files;
  - quality trimming;
  - heterozygous SNV calling from double peaks;
  - heterozygous indel detection by allele separation;
  - automatic handling of reverse reads;
  - per-variant chromatograms;
  - upload in the Add-patient pop-up.
- **All DNA queries:** reverse-orientation sequences are detected and reverse-complemented automatically.
- **Benchmark tool** (`benchmark` command) against ClinGen Evidence Repository expert-panel classifications: leave-one-out, ClinVar override off, agreement/critical discordance/confusion matrix/per-criterion precision-recall.
- **First-run setup** now also installs Biopython for faster alignment.
- **Fixed:** exon number for cDNA references was always reported as exon 1.
- **Examples:** `examples/family/family_example.csv`, `examples/traces/` (simulated .ab1 files).
- **Tests:** 29 automated tests (family, traces, benchmark added).
- Version 8.0.0.

## V7 (2026-11-19)
- **Known variants are now classified from ClinVar at any review level** by default (previously ≥2★ only).
  - For example, SOD1 p.Cys7Ser (ClinVar Pathogenic, 1★) is now reported as Pathogenic instead of VUS.
  - 0-1★ records are labelled "limited review". The automated ACMG class and criteria are always shown next to them.
- **New setting "Classify known variants from ClinVar"** (any / ≥1★ / ≥2★ / ≥3★ / off); `--clinvar-min-stars`, `--no-clinvar-class` on the command line (`run` and `cohort`).
- **Your own known-variant list** is also used as a classification source and takes priority over ClinVar for the same variant.
- **Conflicting ClinVar classifications** are not used; the automated class is reported with a note to review the submissions.
- **Protein changes that two nucleotide changes could explain** (e.g. Cys7Ser: c.19T>A or c.20G>C) are matched to ClinVar through all candidates.
- **SOD1 example:** Patient 4 (p.Cys7Ser) added; new regression tests.
- Version 7.0.0.

## V6 (2026-11-19)
- **Fixed:** known ClinVar pathogenic variants given as protein sequences (e.g. SOD1 p.Lys4Glu, p.Ala5Ser) were reported as VUS.
  - **ClinVar review stars** are now read. A reviewed record (≥2★, no conflicts) gives the **reported class**. The automated ACMG class, criteria and basis are always shown next to it. This can be switched off in the window or with `--no-clinvar-class`.
  - **Protein input:** the nucleotide change is inferred from the RefSeq codon (e.g. p.Lys4Glu → c.10A>G). It is used for gnomAD (PM2), REVEL (PP3) and exact ClinVar matching.
  - **PS1 for protein input** applies when the same amino-acid change is an established pathogenic variant. PM5 is no longer combined with PS1.
- **Fixed:** online search "not available":
  - **NP_ protein accessions** in the reference header (e.g. `NP_000445.1`) are resolved to the gene and transcript via NCBI.
  - **SSL certificates** now work with python.org macOS builds (certifi, or the system bundle). The start files install `certifi`, including into environments set up by earlier versions.
  - **Clear error messages** explain SSL, firewall/proxy, no internet and timeouts.
  - **New "Check connection" button** and `check-connection` command.
- **Fixed:** the review status "...no conflicts" was misread as conflicting.
- **Inheritance from UniProt** now records both modes when a disease is AD and AR.
- **SOD1** added to the built-in knowledge (ALS1, STAHP).
- **Clearer variant sentences** for protein input, which previously read "Not missense...".
- **Report layout fix:** the criteria table no longer stretches.
- **SOD1 protein example** added (`examples/SOD1_protein/`); regression tests now cover this case.
- Version 6.0.0.

## V5 (2026-11-19)
- **Multiple patients in one analysis:**
  - "+ Add patient" opens a pop-up window for ID, sex, age, indication, phenotypes, reviewer evidence, the patient's own sequence(s) and notes;
  - "Save & add another" stays open for the next patient;
  - each patient can be edited, copied or removed;
  - duplicate IDs and missing sequences are rejected.
- **Import a patient list** from CSV/TSV/JSON (`examples/patients_example.csv`).
- **Results:**
  - new patient selector;
  - new "Cohort overview" tab showing every patient and variant.
- **Run folder layout:**
  - one sub-folder per patient;
  - `cohort_summary.html` / `.tsv`, `run_summary.json`, `shared_inputs/`.
- **One patient's failed analysis no longer stops the others;** the reason is reported.
- **Faster online lookups:** they are shared across patients in a run.
- **History** lists all patients of each run. Single-patient runs from V3/V4 still open.
- **Command line:** new `cohort` command (`--patients table.csv`).
- **Example:** the HBB example now loads 5 patients, each with their own clinical details.
- Version 5.0.0.

## V4 (2026-11-19)
- **Choose where results are saved:**
  - "Choose folder..." opens a native folder dialog (Finder on macOS, folder browser on Windows); you can also type a path.
  - The choice is remembered.
  - The default is `Documents/ACMG-VarCall Results`.
- **Automatic annotation (default on), built from public databases:**
  - ClinVar: known P/LP variants, hotspots, LOF mechanism;
  - NCBI RefSeq: transcript sequence, CDS and exons;
  - UniProt: domains, functional sites, diseases and inheritance;
  - gnomAD: constraint;
  - Ensembl VEP, for each detected variant: gnomAD allele frequency, REVEL, SpliceAI, AlphaMissense and rsID.
  - It degrades gracefully when offline.
- **References can be given as a RefSeq accession** (e.g. `NM_000518.5`); the sequence is downloaded.
- **Known pathogenic variants:** the automatic ClinVar search plus a manual list or file (window box, `-k` on the command line).
- **New "Preview online annotation" button.**
- **"Annotation sources" section** in every report, recording what was retrieved and when.
- **Responses are cached for 30 days** in `~/.acmg_varcall/cache`, for reproducibility.
- **PM1** now also uses UniProt active/binding sites and ClinVar-derived hotspots.
- **Command line:** new `--offline` and `-k/--known` options.
- Version 4.0.0.

## V3 (2026-11-19)
- **Opens as its own application window instead of a Safari or browser tab:**
  - a native window via pywebview (installed automatically on first start);
  - otherwise an Edge/Chrome app-mode window;
  - the browser is used only as a last resort.
- **macOS app** (`ACMG-VarCall.app`, with an icon) that starts without a Terminal window.
- **Windows:** the start file runs without leaving a console window open. New `Create_Desktop_Shortcut_Windows.bat` adds a Desktop icon.
- **New Results screen showing all outputs inside the window:**
  - variant maps, full report, alignments, clinical summary and file list;
  - buttons to open the results folder and to print or save a PDF.
- **History screen** to re-open previous analyses.
- **PNG images** of every variant map are now saved next to the SVG files. The command line also does this when `cairosvg` or `rsvg-convert` is available.
- New `app` command, with `--window native|app|browser`.
- Version 3.0.0.

## V2 (2026-11-19)
- **Double-click start files:** `Start_Windows.bat`, `Start_Mac.command`, `start_linux.sh`.
- **New browser-based window** (`acmg-varcall web`):
  - paste or load sequences, annotation and patient details;
  - "Load HBB example" button;
  - results table with links to the report, maps, alignments, summaries, TSV and VCF.
  - It uses only the Python standard library and runs on 127.0.0.1, so no data leaves the computer.
- **Each run is saved** to `results/<date_time>_<patientID>/`, together with a copy of its inputs.
- **`acmg-varcall gui`** (tkinter window) now falls back to the browser window if tkinter is missing.
- Version 2.0.0.

## V1 (2026-11-19)
- **Initial release:**
  - sequence input;
  - local alignment and variant calling;
  - HGVS nomenclature and consequence prediction;
  - ACMG/AMP + ClinGen SVI classification;
  - SVG variant maps;
  - HTML/Markdown/JSON clinical summaries, TSV and VCF;
  - CLI, wizard and tkinter GUI.
