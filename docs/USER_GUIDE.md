# CLAVAGEN - user guide

Clinical variant calling from **DNA or protein sequences**. It compares a patient sequence (Sanger consensus, amplicon, assembled read, cDNA or a protein sequence) with a reference, calls and normalises the variants, and annotates them with HGVS nomenclature and their molecular consequence. Each variant is then classified with the **ACMG/AMP 2015 guidelines** plus the ClinGen SVI refinements. For every sequence you get a variant map, a local alignment, an interactive HTML report and a clinical summary.

> **Decision-support tool only.** All classifications must be reviewed by a qualified clinical scientist or medical geneticist before clinical use. Some criteria need external data (de novo status, segregation, functional studies, case-control data), and the software accepts these as reviewer evidence.

Pure Python 3.8+, with **no required dependencies**. Install `biopython` for faster alignment of long sequences.

## Start the program (no installation needed)

The only requirement is **Python 3.8 or newer**, from https://www.python.org/downloads/. On Windows, tick **"Add Python to PATH"** during installation.

| System | How to start | First time only |
|---|---|---|
| **macOS** | Double-click **`CLAVAGEN.app`** (or `Start_Mac.command`). Keep the app inside this folder; you can drag an alias of it to the Dock. | Right-click → **Open** → **Open**, to pass the Gatekeeper check for unsigned apps. |
| **Windows** | Double-click **`Start_Windows.bat`**. Run `Create_Desktop_Shortcut_Windows.bat` once to get a Desktop icon. | If SmartScreen warns you: *More info → Run anyway*. |
| **Linux** | Run `./start_linux.sh`. | - |

On the first start, a private environment is set up inside this folder (`.venv`). It installs the native-window component (`pywebview`), which takes about a minute and needs internet once. After that the program starts immediately.

**The program opens in its own application window, not in Safari or another browser tab.** It tries these in order:
1. a native window: macOS WebKit, Windows WebView2, or GTK/Qt on Linux;
2. an Edge/Chrome/Brave "app mode" window, which has no address bar;
3. only as a last resort, the default browser.

The interface is served locally on `127.0.0.1`; no data leaves the computer. Closing the window quits the program. If something goes wrong on a Mac, there is a log at `~/Library/Logs/CLAVAGEN.log`.

### Several genes in one run (reference list and library)
Section **1. Reference genes** holds any number of references, one per gene. You can add them:
- from **FASTA files** (several at once, multi-FASTA allowed);
- as a **RefSeq accession** (e.g. `NM_000518.5`, downloaded at run time);
- by **pasting** a sequence;
- from your **reference library**.

Each row has its own settings:

| Setting | Meaning |
|---|---|
| Gene | Detected from the header (`(HBB)`, `NM_` accession); editable |
| Inheritance | *from databases* (default), AD, AR, AD/AR, XLR, XLD. For example, set SOD1 to **AD/AR** to analyse the recessive p.Asp91Ala family correctly. |
| CDS range | 1-based coding region if it can't be detected automatically |
| Preview | Online annotation for that gene (*Preview online annotation* previews all of them) |

**How patients are matched:** each patient's sequence is matched to the best-fitting reference automatically. To force a gene, use *Reference / gene* in the *Add patient* pop-up, or a `reference` column in the patient CSV.

**Library:** **Save to library** stores the current references, with their settings, for future sessions. Add them back with **Add from library**. They are kept in `~/.clavagen/reference_library.json`.

**Example:** **Load multi-gene example** loads HBB and SOD1 with two families (12 people).

**Command line:** give several `-r` options and set per-reference options with `--ref-setting`, keyed by gene symbol (`SOD1:...`) or reference ID (`NM_000454.5:...`):
```
python3 -m clavagen cohort --patients patients.csv -r HBB.fa -r SOD1.fa --ref-setting SOD1:inheritance=AD/AR -o results
```
Settings are matched to references by ID or gene symbol, and `NM_000454` also matches `NM_000454.5`. In an annotation file they go under `"references": {"<ID>": {...}}`.

### Sanger trace files (.ab1)
Patients can be analysed directly from Sanger chromatograms. Attach `.ab1` files in the *Add patient* pop-up (forward and/or reverse reads), or give them on the command line (`-q read.ab1`). Several reads per patient are allowed. The software:
- **Reads the file:** it reads the ABIF format directly; this was checked against Biopython on three real instrument files.
- **Trims low-quality ends:** it uses Mott trimming on the basecaller quality values. If a file has no quality values, quality is estimated from peak purity.
- **Calls heterozygous SNVs from double peaks:** a secondary peak at least 30% of the main peak counts, but only if it is a genuine peak and not the tail of a neighbouring one.
- **Detects heterozygous insertions and deletions:** after one, the two alleles are out of register and the trace shows sustained double peaks. The trace is split into the two alleles, and the second allele is aligned to the reference. For example, it reports a heterozygous `c.126_129del`.
- **Handles reverse reads automatically:** they are reverse-complemented.
- **Shows a chromatogram for each variant** in the report, and saves it as SVG and PNG.

Simulated example traces are in `examples/traces/`.

### Families and couples
Give relatives the same **Family ID**, together with their **relationship to the proband** (proband, mother, father, sibling, child, partner, other) and **affected status**. You can do this in the *Add patient* pop-up or with the CSV columns `family, relationship, affected, parentage_confirmed`. The family is then analysed together:

| Situation | What the software does |
|---|---|
| Variant absent in both parents | **De novo**: PS2 if maternity and paternity are confirmed, otherwise PM6 |
| Two variants in the same gene in the proband | **Phase** from the parents' genotypes. In trans with a P/LP variant in a recessive gene gives PM3 (Moderate if both parents were tested, Supporting if one); in cis with a P/LP variant gives BP2. Compound heterozygosity is reported. |
| Affected relatives carry the variant | **PP1** (simplified count: 2–3 Supporting, 4–5 Moderate, ≥6 Strong; check the informative meioses) |
| An affected relative does not carry it (dominant gene) | **BS4** (consider phenocopy) |
| P/LP variant inherited from an unaffected parent (dominant gene) | Reduced-penetrance note |
| Proband + partner, mother + father | **Couple reproductive risk** (e.g. HbS × HbC → 25% risk of an affected child; dominant → 50%; X-linked → 50% of sons) |

Family evidence is applied to the variant in every family member who carries it, and the members are re-classified. Your own reviewer evidence still takes priority.

The results include a **pedigree drawing** with genotypes, a **segregation table**, the findings, and the couple risks. These appear in the cohort overview and in `cohort_summary.html`, with `family_<ID>_pedigree.svg/png` saved alongside. Each family member's report also gets a *Family analysis* section. To try it, use **Load family example** in the window, or `examples/family/family_example.csv`.

### Benchmark against ClinGen expert panels
```
python3 -m clavagen benchmark -o benchmark/            # downloads the ClinGen Evidence Repository (internet needed)
python3 -m clavagen benchmark --erepo erepo.tsv --max 300 --genes BRCA1,MYH7,PAH
```
**How it works:**
- Variants come from the ClinGen Evidence Repository (expert-panel classifications); each needs an exonic `NM_` c. description.
- For each one, a patient sequence carrying the variant is built from the RefSeq transcript (±150 nt) and analysed exactly like user input.
- The **ClinVar-based reported class is switched off**, and the variant's own known-variant record is removed ("leave-one-out"), so the benchmark measures the automated ACMG engine only.
- **Gene-level evidence is recounted without the tested variant.** Hotspots (residues with ≥2 P/LP missense in ClinVar) and the loss-of-function mechanism (≥3 P/LP truncating variants) are derived from ClinVar, so the tested variant could otherwise count towards its own PM1 or PVS1. Values from the built-in gene list, your annotation file or gnomAD constraint are kept.

**Outputs:**
- `benchmark_report.html`, which reports:
  - 3-tier and exact 5-class agreement;
  - **critical discordance** (P/LP vs B/LB);
  - a confusion matrix;
  - per-criterion TP/FP/FN, precision and recall against the experts' applied criteria;
  - agreement by gene;
  - **why each classification differs**, with every discordant variant given one main cause:
    - *missing evidence:* the expert classification, recalculated from only the criteria a sequence-based tool can compute, equals the tool's class. The experts used functional, case, segregation or phenotype data.
    - *criterion assignment:* the tool applied a computable criterion the experts did not, or missed one they applied.
    - *evidence strength:* same computable criteria, different strength (e.g. PVS1_Strong vs PVS1).
  - agreement with the **expert class recalculated from computable criteria only**, which measures errors of the tool rather than missing information.
- `benchmark_results.tsv`, with one row per variant, including the recalculated expert class, the cause of any disagreement, and the criteria applied only by the tool, only by the experts, or with a different strength.
- `benchmark_summary.json`.

**Run time:** each variant needs online look-ups. A 200-variant run takes roughly 15–30 min; Biopython, installed by the start files, keeps the alignment fast.

### Several patients in one analysis
The reference sequence, annotation and known variants are entered once and shared by all patients. Each patient then has their own sequences, details and evidence.
- **Add a patient:** in section **4. Patients & clinical evidence**, click **+ Add patient**. A pop-up window asks for:
  - Patient ID (required, must be unique);
  - sex, age and indication;
  - phenotypes and reviewer evidence;
  - the patient's own sequence(s): paste them or load a file (required);
  - notes.
- **Save & add another** keeps the pop-up open for the next patient.
- **Manage the list:** each patient has **Edit**, **Copy** (useful for family members) and **Remove** buttons.
- **Import patient list(s)** accepts several files at once and can be used repeatedly to add more patients. There is no limit on the number of patients. If an imported list contains **patient IDs that are already in the list** (for example two families that both use `P1_Father`), a dialog asks whether to **rename the new patients** (the family ID is added in front, e.g. `FAM-HBB-01_P1_Father`), **replace** the existing ones, **skip** the duplicates or cancel. If a **Family ID** is already used by another family, it is renamed too (e.g. `FAM1_2`), so different families are never merged by mistake.
- **Import patient list** reads a CSV, TSV or JSON table with the columns `id, sex, age, indication, phenotypes, evidence, query`. `query` is the sequence itself or the name of a FASTA file next to the table. See `examples/patients_example.csv`.
- **Run analysis** analyses every patient. Phenotypes and reviewer evidence apply only to their own patient.
- **Results:**
  - **Cohort overview** shows all patients together.
  - The **Patient** list at the top switches between patients' maps, reports, alignments and files.
- **Folder structure:** each run is one folder:
  - `cohort_summary.html` and `cohort_summary.tsv`, covering all patients and variants;
  - `run_summary.json`;
  - `shared_inputs/`;
  - one sub-folder per patient with the full outputs: `report.html`, maps (SVG/PNG), alignments, `clinical_summary.md`/`.json`, `variants.tsv`/`.vcf` and that patient's inputs.
- **Command line:** `python3 -m clavagen cohort --patients patients.csv -r reference.fa -a annotation.json -o results/`

### Reported class: ClinVar and automated ACMG
- **Automated ACMG class:** each variant gets an automated ACMG/AMP class from the evidence the software can find itself:
  - population frequency, in-silico scores, protein domains and hotspots;
  - other known variants at the same residue;
  - whether loss of function is the disease mechanism.
- **Why it can be lower than ClinVar:** literature evidence, such as functional studies (PS3), case-control data (PS4) and segregation (PP1), is not available to the software. So a well-known pathogenic missense variant can come out as VUS or LP from automation alone.
- **Known variants get the ClinVar class (default):** if a variant has a ClinVar record, its ClinVar classification is **reported**, whatever the review level. A matching entry in your own known-variant list is used the same way, and takes priority.
- **Review stars always shown:** they tell you how much weight the record carries:

  | Stars | Meaning |
  |---|---|
  | 0★ | No assertion criteria; weak evidence |
  | 1★ | Criteria provided, single submitter. Labelled *limited review*. |
  | 2★ | Multiple submitters, no conflicts |
  | 3★ | Expert panel |
  | 4★ | Practice guideline |

- **Always shown alongside:** the automated ACMG class, its criteria, the VCV ID, the condition and the basis.
- **Conflicting records:** ClinVar records with *conflicting classifications* are never used; the automated ACMG class is reported, with a note.
- **Raising the minimum:** use the setting *Classify known variants from ClinVar* in section 2: any record (default), ≥1★, ≥2★, ≥3★, or off (automated ACMG only). From the command line, use `--clinvar-min-stars N` or `--no-clinvar-class`.

### Protein sequences
For a protein query, the software looks up the RefSeq transcript (from an `NM_` or `NP_` accession, or the gene symbol). It then infers the single-nucleotide change that explains each amino-acid change, for example `p.Lys4Glu` → `c.10A>G` (codon AAG→GAG). That inferred change is used to look up gnomAD frequency, REVEL, SpliceAI and the exact ClinVar record. If several nucleotide changes are possible, all of them are checked. Inferred changes are labelled and must be confirmed by DNA sequencing. Zygosity cannot be determined from protein. A worked example is in `examples/SOD1_protein/`.

**Truncating changes in protein queries (V11).** A frameshift is called when the query matches the reference up to a point and then continues with unrelated residues, e.g. `p.(Phe42LeufsTer19)`; the small insertion/deletion that explains the new residues is inferred from the transcript (e.g. `c.126_129del`), or all candidates are listed when codon redundancy allows several. A query that ends in `*` early is called as a premature stop, and residues past the reference stop as a stop-loss. **End a complete protein sequence with `*`**: without it, a query that just stops early is treated as a partial sequence and gets a coverage warning.

### Genotype per patient and compound heterozygosity
All of a patient's files (amplicons, forward and reverse traces, consensus) are combined per gene, so a genotype split over several files is interpreted as a whole. Each variant is labelled by type: null (nonsense, frameshift, canonical splice, start-loss), missense, in-frame, stop-loss, splice region or synonymous.

The genotype is then stated whatever the classification. It is one of:
- heterozygous;
- homozygous;
- **compound heterozygous, in trans** (phase from the parents);
- in cis;
- **phase unknown (possible compound heterozygote)**.

The class of each allele appears next to it, e.g. *"Compound heterozygous, in trans (confirmed by testing both parents) - both variants VUS [null (nonsense) / missense]"*. A comment describes the allele combination (null/null, null/missense, missense/missense).

**Indels in a text or consensus sequence:** these can show only one allele, so their zygosity is reported as *unresolved* until a Sanger trace is used.

### Truncating changes in protein sequences
A protein query is checked for frameshifts, premature stops and stop-loss: e.g. `p.(Phe42LeufsTer19)`, `p.(Gln40Ter)`, or an extension.

**End a complete protein with `*`.** Without it, a query that just stops early is treated as a partial sequence. When a DNA reference or RefSeq transcript is available, the DNA change behind a frameshift is inferred, e.g. `c.126_129del`. If several changes give the same protein, all of them are listed.

### Troubleshooting the online search
Click **Check connection** (next to *Preview online annotation*), or run `python3 -m clavagen check-connection`. It tests NCBI, UniProt, gnomAD and Ensembl and explains each failure:
- **"secure connection failed (SSL certificate check)"** is common on macOS with the python.org Python. Start the program with the start file (V6 installs certificates automatically), or open *Applications → Python 3.x → Install Certificates.command* once.
- **"blocked by a firewall/proxy"**: a hospital or university network may block these sites. Ask IT to allow `eutils.ncbi.nlm.nih.gov`, `rest.uniprot.org`, `gnomad.broadinstitute.org` and `rest.ensembl.org`, or use another network.
- **"no internet connection"**: check the network.
- **"Gene could not be identified"**: type the gene symbol in *Gene symbol*, or use a reference header with an `NM_`/`NP_` accession or `(GENE)`.

### Where results are saved
Results go to **`Documents/CLAVAGEN Results`** by default. To change this, use **Choose folder...** at the bottom of *New analysis*: a standard Finder or Windows folder dialog opens. You can also type a path and click **Use this folder**. The choice is remembered for future sessions. From the command line, use `-o FOLDER`.

### Automatic annotation and known pathogenic variants
This is on by default (the *Search online databases automatically* checkbox; `--offline` on the command line turns it off). The software identifies the gene from the **Gene symbol** box, an `NM_` accession, or `(GENE)` in the reference header. It then builds the annotation itself:

| Source | What is used | ACMG use |
|---|---|---|
| **ClinVar** (NCBI E-utilities) | Pathogenic / likely pathogenic variants in the gene | PS1, PM5, database assertions, hotspots (≥2 P/LP missense at a residue), LOF mechanism (≥3 P/LP truncating variants) |
| **NCBI RefSeq** | Transcript sequence, CDS and exon structure. You can type just `NM_000518.5` as the reference. | Correct c. numbering, NMD prediction (PVS1 strength) |
| **UniProt** | Domains/regions, active and binding sites, function, diseases (MIM, inheritance) | PM1, disease summary, inheritance |
| **gnomAD** | Gene constraint (LOEUF, missense Z) | LOF mechanism support, PP2 (missense Z ≥ 3.09) |
| **Ensembl VEP** | For each detected variant: gnomAD allele frequency, REVEL, SpliceAI, AlphaMissense, dbSNP rsID, ClinVar significance | BA1, BS1, PM2, PP3, BP4, BP7 |

- **Preview online annotation** shows what was found before you run.
- **Order of precedence:** online data < built-in curated knowledge < annotation JSON file < manual known variants.
- **Manual known variants:** add your own (in-house, published or lab database) in the *Add known pathogenic variants manually* box, or load a `.txt`, `.csv`, `.tsv` or `.json` file. Write one variant per line, for example `c.20A>T | p.Glu7Val | Pathogenic | In-house`. A c. or p. change on its own is enough. From the command line, use `-k file.tsv` or `-k "c.20A>T|p.Glu7Val|Pathogenic|source"`.
- **What is sent online:** only gene symbols, transcript accessions and variant names. Patient identifiers and sequences are never sent.
- **Caching:** responses are cached for 30 days in `~/.clavagen/cache`, so a re-analysis is reproducible and also works offline.
- **Record of sources:** every report has an **Annotation sources** section listing what was retrieved, when, and anything that was unavailable.
- **Without internet:** the analysis still runs on the sequence alone (plus any file or manual entries), and the report says which sources were missing.
- **REVEL licence:** REVEL scores are free for non-commercial use only.

### Using the window
- **New analysis:** paste or load the query sequence(s), the reference(s) and the annotation JSON, then add patient details, phenotypes and reviewer evidence. **Load HBB example** fills in the demo.
- **Results** (opens automatically after a run):
  - **Variant maps:** the map images and variant tables for each sequence.
  - **Full report:** ACMG criteria evidence, interpretation and recommendations.
  - **Alignments:** the local alignment for each sequence.
  - **Clinical summary.**
  - **Files:** click a file to open it in its normal program. **Open results folder** and **Print / save PDF** are also here.
- **History:** re-open any earlier run.
- **Saved files:** every run is stored in `<results folder>/<date_time>_<patientID>/`. It contains the maps (SVG and PNG), alignments, `report.html`, `clinical_summary.md`/`.json`, `variants.tsv`/`.vcf` and a copy of the inputs.

Command-line equivalents are `python3 -m clavagen app` (the window) and `python3 -m clavagen run ...` (batch).

## Install (optional, for command-line use)

```bash
pip install .            # or: pip install .[fast]   (adds Biopython)
clavagen --version   # or: python -m clavagen --version
```

## Quick start

```bash
clavagen example -o demo        # writes HBB demo inputs and runs them -> demo/results/report.html

clavagen run \
  -q patients.fa -r HBB_reference.fa -a HBB_annotation.json \
  -p "microcytic anaemia" -e "c.20A>T=PS3,PS4" \
  --patient-id P001 --sex female -o results/
```

Other entry points:

| Command | Purpose |
|---|---|
| `clavagen run` | Full analysis (call → annotate → classify → map → report) |
| `clavagen align` | Local or global alignment only (BLAST-style text output) |
| `clavagen translate` | 3-frame translation, ORF finder, GC content (`--reverse` for reverse complement) |
| `clavagen wizard` | Interactive prompts: paste sequences, no command-line flags needed |
| `clavagen app` | Application window (what the start files open). `--window native/app/browser` chooses the window type. |
| `clavagen web` | The same interface in a normal browser tab |
| `clavagen benchmark` | Benchmark against ClinGen Evidence Repository expert classifications |
| `clavagen check-connection` | Test access to the online databases |

## Inputs

- **Query** (`-q`, repeatable) and **reference** (`-r`, repeatable) can each be a FASTA or multi-FASTA file, `NAME=SEQUENCE`, or a raw sequence.
- The molecule type is detected automatically. Supported combinations:
  - DNA query vs DNA reference (cDNA or genomic with exons).
  - Protein query vs protein reference.
  - Protein query vs DNA reference: the reference CDS is translated first.
- **Heterozygous calls:** two-base IUPAC codes in a DNA query (e.g. `W` = A/T, `R` = A/G) are called as heterozygous SNVs. `N`/`X` are reported as no-calls.
- **Several references:** each query is matched to the best-scoring reference automatically, or you can pair them explicitly with `--pair QUERY_ID:REF_ID`.
- **CDS:** taken from the annotation or `--cds 51-494`. If neither is given, the software uses the whole reference (when it runs ATG…stop in frame) or the longest ORF.

### Annotation (`-a annotation.json`)
See `examples/annotation_template.json`. The annotation can hold:
- Gene and transcript IDs, CDS and exons, or `exon_boundaries_c` for cDNA references (used for NMD prediction).
- Inheritance and LOF mechanism.
- Protein domains, with `critical` flags (used for PM1).
- Hotspots.
- Known variants (PS1, PM5 and database assertions).
- Population allele frequencies (BA1, BS1, PM2).
- REVEL and SpliceAI scores (PP3, BP4, BP7).
- BA1 exceptions and threshold overrides.
- Disease summaries and management notes.

You can also nest any of these under `references.{REF_ID}` to give per-reference settings.

If you don't supply an annotation, `--gene SYMBOL` loads a small built-in starter knowledge base (HBB, CFTR, BRCA1/2, TP53, MLH1, FBN1, LDLR, PAH, DMD, MYH7). Genes on the **ACMG SF v3.3** secondary-findings list (84 genes) are flagged.

### Reviewer evidence (`-e`, `--evidence-file`)
- Add a criterion: `"c.20A>T=PS3,PP1_Strong"`.
- Remove an automated criterion: `"p.Glu7Val=-PM2"`.
- Apply to all variants: `"*=PP4"`.

Strength modifiers are `_VeryStrong`, `_Strong`, `_Moderate` and `_Supporting`. Reviewer evidence always overrides the automated call.

## What is automated

| Criterion | How |
|---|---|
| PVS1 | Nonsense, frameshift, canonical ±1/2 splice, start-loss or exon-boundary change in a gene with `lof_mechanism: true`. Strength follows the ClinGen PVS1 tree: NMD prediction (>50 nt upstream of the last junction), fraction of protein lost, critical region, start-loss → Moderate. |
| PS1 / PM5 | Same amino-acid change, or a different pathogenic missense at the same residue, in `known_variants`. |
| PM1 | Annotated hotspot or critical domain. The PP3+PM1 total is capped at Strong. |
| PM2 / BS1 / BA1 | Allele frequency against inheritance-specific thresholds (configurable). PM2 is applied at Supporting strength (SVI 2020), and BA1 exceptions are respected. |
| PM4 / BP3 | In-frame indel or stop-loss; repeat context is detected for BP3. |
| PP2 / BP1 | Gene-level flags `missense_constrained` and `truncating_only`. |
| PP3 / BP4 | REVEL scores with calibrated thresholds (Pejaver 2022) and SpliceAI (≥0.2 / ≤0.1). If no score is supplied, a transparent **heuristic** is used instead: BLOSUM62, Grantham, physico-chemical class and ortholog conservation. It is flagged, and `--strict` disables it. |
| BP7 | Synonymous or deep-intronic variant outside the splice region. |
| PP4 | Patient phenotype matched against the gene's disease spectrum (heuristic, flagged). |
| PS2–PS4, PM3, PM6, PP1, BS2–BS4, BP2, BP5 | Reviewer evidence only. A PM3 reminder is added for recessive (AR) genes. |

**Combining the criteria:**
- The ACMG/AMP Table 5 rules decide the classification, with one ClinGen SVI addition: PVS1 + a Supporting criterion → Likely pathogenic.
- The Tavtigian **points** total is reported alongside. VUS results get a sub-tier (VUS-high / mid / low).
- When the rules and the points disagree, the report flags it for the reviewer.
- PP5 and BP6 are not used, following ClinGen SVI guidance.

## Outputs (in `-o` folder)

| File | Content |
|---|---|
| `<query>__vs__<ref>_map.svg` / `.png` | Variant map for each sequence, with:<br>• lollipops coloured by ACMG class, stem and ring coloured by variant type<br>• CDS, UTR, exon junctions and domains<br>• query coverage and no-calls<br>• a zoom panel per variant: ref vs query bases, codon change, BLOSUM62/Grantham |
| `<query>__vs__<ref>_alignment.txt` | Local alignment in BLAST style (score, identity, similarity, gaps, spans) |
| `report.html` | Self-contained interactive report:<br>• summary table filterable by class<br>• clinical interpretation, condition summaries and recommendations<br>• the maps and colour-coded alignments<br>• criterion-by-criterion evidence table<br>• QC, traceability (input SHA-256 hashes) and disclaimer |
| `clinical_summary.md` / `.json` | Clinical features, scores, criteria, a short disease summary, genotype interpretation and recommendations |
| `variants.tsv` | One row per variant (HGVS c./p., consequence, metrics, class, points, criteria) |
| `variants.vcf` | VCF 4.2 with GT and ACMG INFO fields. Use `chrom` and `genomic_offset` in the annotation to give genomic coordinates. |

## Features for clinicians

- **Genotype-aware interpretation:**
  - AR genes: carrier status, biallelic results, or possible compound heterozygosity (with a phasing reminder).
  - X-linked genes: interpretation by sex.
  - AD genes: diagnosis or predisposition.
- **Secondary findings:** genes on the ACMG SF v3.3 list are flagged.
- **QC:**
  - Warnings for low identity, partial CDS coverage, no-calls, unresolved ambiguity codes and high variant density.
  - CDS sanity checks (ATG start, in-frame, internal stops).
- **Nomenclature:** HGVS 3′-normalisation of indels, duplication detection, and intronic c. numbering (e.g. `c.92+1G>A`).
- **Amino-acid metrics:** BLOSUM62, Grantham distance and class, hydropathy change and property class change.
- **Conservation:** supply ortholog proteins with `--orthologs` to get per-residue conservation.
- **Traceability:** each report records the software version, run date, settings and input hashes.

## Tests

```bash
pip install pytest && pytest -q
```

## Limitations

- The software analyses only the region covered by the query. CNVs, large rearrangements and mosaicism are out of scope.
- Zygosity comes from IUPAC codes. A consensus sequence cannot show heterozygous indels.
- The built-in knowledge base and the example population and in-silico values are illustrative. Use current gnomAD, ClinVar, ClinGen and OMIM data in your annotation.
- Disease-specific ClinGen VCEP specifications override the generic thresholds: set them with `thresholds` and reviewer evidence.
