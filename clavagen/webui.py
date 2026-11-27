"""CLAVAGEN desktop application window.

The interface is served by a small local server (127.0.0.1 only - no data leaves
the computer) and displayed in its own application window:

1. a native window via ``pywebview`` (installed automatically by the start files), or
2. a chrome-less "app mode" window of Edge / Chrome / Chromium / Brave, or
3. as a last resort, a tab in the default web browser.
"""

from __future__ import annotations

import base64
import datetime as _dt
import json
import mimetypes
import os
import re
import shutil
import socket
import subprocess
import sys
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, unquote, urlparse

from . import __version__

APP_TITLE = "CLAVAGEN"

PAGE = r"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>CLAVAGEN</title>
<style>
*{box-sizing:border-box}
html,body{height:100%;margin:0}
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;background:#f4f6f9;color:#1d2433;display:flex;flex-direction:column}
.bar{background:#13294b;color:#fff;display:flex;align-items:center;padding:0 18px;height:52px;flex:none}
.bar .logo{font-weight:700;font-size:17px;margin-right:28px}.bar .logo small{font-weight:400;opacity:.7;font-size:12px;margin-left:6px}
.nav button{background:transparent;color:#cfd8e6;border:0;padding:16px 14px;font-size:14px;cursor:pointer;border-bottom:3px solid transparent}
.nav button.on{color:#fff;border-bottom-color:#f0b400}.nav button:disabled{opacity:.4;cursor:default}
.bar .right{margin-left:auto;font-size:12px;opacity:.7}
.view{display:none;flex:1;overflow:auto}.view.on{display:block}
#v-results.on{display:flex;flex-direction:column;overflow:hidden}
.wrap{max-width:1180px;margin:0 auto;padding:18px 22px 40px}
.card{background:#fff;border:1px solid #e3e7ee;border-radius:10px;padding:16px 20px;margin-bottom:16px}
h2{font-size:15px;margin:0 0 10px}
label{font-size:13px;font-weight:600;display:block;margin:10px 0 4px}.hint{font-weight:400;color:#6b7280;font-size:12px}
textarea,input[type=text],input[type=number],select{width:100%;border:1px solid #cfd6e0;border-radius:6px;padding:7px 9px;font-size:13px;background:#fff}
textarea{font-family:Menlo,Consolas,monospace;min-height:120px}
.row{display:grid;grid-template-columns:1fr 1fr;gap:16px}.row3{display:grid;grid-template-columns:1fr 1fr 1fr;gap:16px}
.btns{display:flex;gap:10px;align-items:center;margin-top:14px;flex-wrap:wrap}
button.b{border:0;border-radius:6px;padding:9px 18px;font-size:14px;cursor:pointer}
.primary{background:#13294b;color:#fff}.secondary{background:#e7ebf2;color:#13294b}.small{padding:6px 12px!important;font-size:13px!important}
.file{font-size:12px;margin-top:4px}
#status{font-size:13px;color:#6b7280}
table{border-collapse:collapse;width:100%;font-size:13px}td,th{border-bottom:1px solid #eee;padding:6px 8px;text-align:left;vertical-align:top}
th{background:#f6f7f9}
.badge{color:#fff;padding:2px 8px;border-radius:10px;font-size:12px;white-space:nowrap}
.warn{color:#b45309;font-size:13px}
.toolbar{background:#fff;border-bottom:1px solid #e3e7ee;padding:10px 18px;flex:none}
.toolbar .top{display:flex;align-items:center;gap:10px;flex-wrap:wrap}
.toolbar .title{font-weight:700;font-size:15px;margin-right:auto}
.chips{display:flex;gap:8px;flex-wrap:wrap;margin-top:8px}
.chip{background:#f4f6f9;border-radius:8px;padding:5px 10px;font-size:12px}
.tabs{display:flex;gap:4px;margin-top:10px}
.tabs button{border:1px solid #d6dce6;background:#f7f8fa;padding:6px 14px;border-radius:6px 6px 0 0;cursor:pointer;font-size:13px}
.tabs button.on{background:#13294b;color:#fff;border-color:#13294b}
.pane{flex:1;overflow:auto;display:none;background:#fff}.pane.on{display:block}
#p-report iframe{width:100%;height:100%;border:0;display:block}#p-report.on{display:flex}
.mapbox{border:1px solid #e3e7ee;border-radius:8px;margin:14px 18px;padding:10px;background:#fff}
.mapbox h3{margin:2px 4px 8px;font-size:14px}.mapbox img{width:100%;height:auto;display:block}
pre.aln{background:#0f172a;color:#e2e8f0;padding:12px;border-radius:8px;font-size:12px;line-height:1.35;margin:14px 18px;overflow:auto}
pre.md{white-space:pre-wrap;font-family:Menlo,Consolas,monospace;font-size:12.5px;margin:14px 18px;line-height:1.45}
.empty{color:#6b7280;padding:60px;text-align:center}
.hist tr{cursor:pointer}.hist tr:hover td{background:#f3f6fb}
.help li{margin:6px 0;font-size:14px}
@media(max-width:800px){.row,.row3{grid-template-columns:1fr}}
.modal-bg{position:fixed;inset:0;background:rgba(15,23,42,.55);display:none;align-items:center;justify-content:center;z-index:50}
.modal-bg.on{display:flex}
.modal{background:#fff;border-radius:12px;width:min(860px,94vw);max-height:92vh;display:flex;flex-direction:column;box-shadow:0 20px 60px rgba(0,0,0,.35)}
.mhead{display:flex;align-items:center;padding:14px 20px;border-bottom:1px solid #e3e7ee;font-weight:700;font-size:16px}
.mhead .x{margin-left:auto;border:0;background:none;font-size:24px;cursor:pointer;color:#6b7280}
.mbody{padding:6px 20px 14px;overflow:auto}.mfoot{display:flex;justify-content:flex-end;gap:10px;padding:12px 20px;border-top:1px solid #e3e7ee}
.ptable td{vertical-align:middle}.lnk{border:0;background:none;color:#1d4ed8;cursor:pointer;font-size:13px;padding:2px 4px}.lnk.red{color:#b91c1c}
.mono{font-family:Menlo,Consolas,monospace}
.famrow{background:#f6f8fb;border-radius:8px;padding:4px 12px 10px;margin-top:10px}.famhead{font-weight:700;font-size:13px;margin-top:6px}
.family{border:1px solid #e3e7ee;border-radius:10px;padding:10px 16px;margin:16px 0;background:#fff}.family img{max-width:100%}
li.risk{color:#b2182b}
</style></head><body>
<div class="bar"><div class="logo">CLAVAGEN<small>v__VERSION__</small></div>
<div class="nav"><button data-v="input" class="on">New analysis</button><button data-v="results" id="nav-results" disabled>Results</button>
<button data-v="history">History</button><button data-v="help">Help</button></div>
<div class="right">Local only &middot; decision support - expert review required</div></div>

<div class="view on" id="v-input"><div class="wrap">
<div class="card"><h2>1. Reference genes <span class="hint" id="rcount"></span></h2>
 <p class="hint" style="margin-top:0">Add one reference per gene (cDNA, genomic or protein FASTA, or a RefSeq accession). Several genes can be
 analysed in one run - each patient's sequence is matched to the right reference automatically.</p>
 <div class="btns" style="margin:0 0 10px">
  <button class="b primary small" id="ref-add-file">+ Add reference file(s)...</button>
  <input type="file" id="ref-file" multiple accept=".fa,.fasta,.fna,.faa,.txt,.seq,.fas,.gb" style="display:none">
  <input type="text" id="ref-acc" placeholder="e.g. NM_000518.5" style="width:160px"><button class="b secondary small" id="ref-add-acc">Add accession</button>
  <button class="b secondary small" id="ref-paste-toggle">Paste sequence...</button>
  <select id="ref-lib" style="width:auto;max-width:230px"><option value="">Add from library...</option></select>
  <button class="b secondary small" id="ref-save-lib" title="Remember these references for future sessions">Save to library</button>
 </div>
 <div id="ref-paste" style="display:none;margin-bottom:10px">
  <textarea id="reference" style="min-height:80px" placeholder=">NM_000518.5 (HBB) ...&#10;ATGGTGCAT...   (one or several FASTA records)"></textarea>
  <button class="b secondary small" id="ref-add-paste" style="margin-top:6px">Add pasted reference(s)</button></div>
 <table class="ptable"><thead><tr><th>Reference ID</th><th>Gene</th><th>Type / length</th><th>Inheritance</th><th>CDS range <span class="hint">(1-based)</span></th><th></th></tr></thead>
 <tbody id="ref-rows"></tbody></table>
 <div id="ref-msg" class="hint" style="margin-top:6px"></div>
</div>
<div class="card"><h2>2. Automatic annotation &amp; known pathogenic variants</h2>
<label style="font-weight:400"><input type="checkbox" id="online" checked> <b>Search online databases automatically</b>
 <span class="hint">- ClinVar (known P/LP variants), NCBI RefSeq (CDS/exons), UniProt (domains, disease), gnomAD (constraint), Ensembl VEP
 (population frequency, REVEL, SpliceAI, AlphaMissense). Only the gene, transcript and variant names are sent - never patient data.</span></label>
<label style="font-weight:400;margin-top:6px"><b>Classify known variants from ClinVar / your known-variant list:</b>
 <select id="clinvar_min" style="width:auto;margin-left:6px">
  <option value="0" selected>Yes - any ClinVar record (all review levels)</option>
  <option value="1">Yes - ClinVar &ge;1&#9733; (criteria provided)</option>
  <option value="2">Yes - ClinVar &ge;2&#9733; (multiple submitters / expert panel)</option>
  <option value="3">Yes - ClinVar &ge;3&#9733; (expert panel / guideline)</option>
  <option value="off">No - report the automated ACMG class only</option></select>
 <span class="hint">- the automated ACMG class, criteria and the ClinVar stars are always shown. Records with conflicting
 classifications keep the automated class.</span></label>
<div class="btns" style="margin-top:6px"><button class="b secondary small" id="preview">Preview online annotation</button>
<button class="b secondary small" id="checkconn">Check connection</button>
<span id="pstatus" class="hint"></span></div>
<div id="preview-box" style="display:none;margin-top:8px;font-size:13px;background:#f6f8fb;border-radius:8px;padding:10px 14px"></div>
<label>Add known pathogenic variants manually <span class="hint">(optional; added to the online list and take priority).
 One per line: <code>c.20A&gt;T | p.Glu7Val | Pathogenic | source</code> - columns separated by | , ; or tab; c. or p. alone is fine</span></label>
<textarea id="known" style="min-height:70px" placeholder="c.20A>T | p.Glu7Val | Pathogenic | In-house database&#10;c.118C>T | p.Gln40Ter | Pathogenic | Published case"></textarea>
<input class="file" type="file" data-target="known" accept=".txt,.csv,.tsv,.json">
</div>
<div class="card"><h2>3. Annotation file <span class="hint">(optional - overrides online values; gene / inheritance / CDS can be set per reference in section 1)</span></h2>
<label>Annotation JSON <span class="hint">gene, CDS/exons, inheritance, domains, known variants, population AF, REVEL/SpliceAI - see examples/annotation_template.json</span></label>
<textarea id="annotation" style="min-height:80px" placeholder='{"gene": "HBB", "cds": [1, 444], ...}'></textarea>
<input class="file" type="file" data-target="annotation" accept=".json">
<label>Gene rule profile <span class="hint">(optional) gene-specific settings from an expert panel specification
 (ClinGen CSpec Registry): thresholds, criteria not used, fixed strengths and caps. Get a template with <code>clavagen evidence-template</code>.
 Copy the values from the published specification - do not guess them.</span></label>
<textarea id="rules" style="min-height:60px" placeholder='{"MYH7": {"source": "ClinGen Cardiomyopathy VCEP v2 (URL)", "thresholds": {"ba1": 0.001}, "disabled": ["PP2"]}}'></textarea>
<input class="file" type="file" data-target="rules" accept=".json">
<div class="row3">
<input type="hidden" id="gene"><input type="hidden" id="cds">
 <div><label>Mode</label><select id="strict"><option value="0">Standard (heuristics flagged)</option><option value="1">Strict (no heuristic evidence)</option></select></div>
</div></div>
<div class="card"><h2>4. Patients &amp; clinical evidence <span class="hint" id="pcount"></span></h2>
 <p class="hint" style="margin-top:0">Each patient has their own sequence(s), details, phenotypes and reviewer evidence. All patients are analysed in one run.</p>
 <div class="btns" style="margin:0 0 10px"><button class="b primary small" id="add-patient">+ Add patient</button>
  <button class="b secondary small" id="import-btn">Import patient list(s) (CSV/TSV)...</button>
  <input type="file" id="import-file" multiple accept=".csv,.tsv,.txt,.json" style="display:none">
  <button class="b secondary small" id="clear-patients">Remove all</button></div>
 <table class="ptable"><thead><tr><th>#</th><th>Patient ID</th><th>Family</th><th>Sex</th><th>Age</th><th>Indication / phenotypes</th><th>Sequence(s)</th><th>Clinical / reviewer evidence</th><th></th></tr></thead>
 <tbody id="patient-rows"></tbody></table>
</div>
<div class="card"><h2>5. Save &amp; run</h2>
<label>Save results in</label>
<div style="display:flex;gap:8px"><input type="text" id="results_dir" style="flex:1"><button class="b secondary small" id="choose">Choose folder...</button>
<button class="b secondary small" id="setdir">Use this folder</button></div>
<div class="btns"><button class="b primary" id="run">Run analysis</button><button class="b secondary" id="example">Load HBB example</button><button class="b secondary" id="famexample">Load family example</button><button class="b secondary" id="multiexample">Load multi-gene example</button>
<button class="b secondary" id="clear">Clear all</button><span id="status"></span></div></div>
</div></div>

<div class="modal-bg" id="modal"><div class="modal">
 <div class="mhead"><span id="m-title">Add patient</span><button class="x" id="m-close" title="Close">&times;</button></div>
 <div class="mbody">
  <div class="row3"><div><label>Patient ID *</label><input type="text" id="m-id" placeholder="e.g. P-2026-014"></div>
   <div><label>Sex</label><select id="m-sex"><option value="">-</option><option>female</option><option>male</option><option>unknown</option></select></div>
   <div><label>Age</label><input type="text" id="m-age" placeholder="e.g. 34"></div></div>
  <div class="famrow"><div class="famhead">Family / couple <span class="hint">(optional - give relatives the same Family ID to analyse them together)</span></div>
   <div class="row3"><div><label>Family ID</label><input type="text" id="m-family" list="fam-list" placeholder="e.g. FAM-014"><datalist id="fam-list"></datalist></div>
    <div><label>Relationship to proband</label><select id="m-relationship"><option value="">-</option><option>proband</option><option>mother</option>
     <option>father</option><option>sibling</option><option>child</option><option>partner</option><option>other</option></select></div>
    <div><label>Affected status</label><select id="m-affected"><option value="unknown">unknown</option><option>affected</option><option>unaffected</option></select></div></div>
   <label style="font-weight:400"><input type="checkbox" id="m-parentage"> Maternity and paternity confirmed <span class="hint">(proband; upgrades de novo evidence from PM6 to PS2)</span></label>
  </div>
  <label>Reference / gene <span class="hint">(leave on auto-detect unless you want to force a gene)</span></label>
  <select id="m-reference"><option value="">Auto-detect (best-matching reference)</option></select>
  <label>Indication</label><input type="text" id="m-indication" placeholder="Reason for testing">
  <label>Phenotypes <span class="hint">comma-separated (used for PP4), e.g. microcytic anaemia, splenomegaly</span></label><input type="text" id="m-phenotypes">
  <details class="famrow" id="m-ce-box"><summary class="famhead" style="cursor:pointer">Clinical evidence (guided) <span class="hint">- describe what you observed;
   CLAVAGEN converts it to ACMG criteria and strengths with the ClinGen rules and shows the calculation in the report</span></summary>
   <label>Applies to variant <span class="hint">(c. or p.; leave empty for all variants found in this patient)</span></label><input type="text" id="ce-key" placeholder="e.g. c.20A>T">
   <div class="row3"><div><label>Informative meioses <span class="hint">(segregation, PP1)</span></label><input type="number" min="0" id="ce-meioses"></div>
    <div><label>Families</label><input type="number" min="0" id="ce-families"></div>
    <div><label>Healthy adults with the genotype <span class="hint">(BS2)</span></label><input type="number" min="0" id="ce-healthy"></div></div>
   <div class="famhead" style="margin-top:8px">Recessive disease: seen in trans with (PM3)</div>
   <div class="row3"><div><label>Pathogenic, phase confirmed</label><input type="number" min="0" id="ce-tp-c"></div>
    <div><label>Pathogenic, phase unknown</label><input type="number" min="0" id="ce-tp-u"></div>
    <div><label>Likely pathogenic, phase confirmed</label><input type="number" min="0" id="ce-tlp-c"></div></div>
   <div class="row3"><div><label>Likely pathogenic, phase unknown</label><input type="number" min="0" id="ce-tlp-u"></div>
    <div><label>Homozygous affected probands</label><input type="number" min="0" id="ce-hom"></div><div></div></div>
   <div class="famhead" style="margin-top:8px">Functional assay (PS3 / BS3)</div>
   <div class="row3"><div><label>Result</label><select id="ce-func"><option value="">none</option><option value="damaging">damaging / abnormal</option><option value="normal">normal / no effect</option></select></div>
    <div><label>OddsPath <span class="hint">(if the assay is validated)</span></label><input type="number" step="any" min="0" id="ce-oddspath"></div>
    <div><label>or strength</label><select id="ce-fstr"><option value="">-</option><option>Supporting</option><option>Moderate</option><option>Strong</option></select></div></div>
   <div class="famhead" style="margin-top:8px">De novo (PS2 / PM6)</div>
   <div class="row3"><div><label>Probands, confirmed de novo</label><input type="number" min="0" id="ce-dn-c"></div>
    <div><label>Probands, assumed de novo</label><input type="number" min="0" id="ce-dn-a"></div>
    <div><label>Phenotype</label><select id="ce-dn-ph"><option value="specific">highly specific for the gene</option><option value="consistent">consistent, not specific</option><option value="heterogeneous">consistent, high heterogeneity</option></select></div></div>
   <div class="famhead" style="margin-top:8px">Case data (PS4)</div>
   <div class="row3"><div><label>Unrelated affected probands</label><input type="number" min="0" id="ce-ap"></div>
    <div><label>Strength per panel rules</label><select id="ce-apstr"><option>Supporting</option><option>Moderate</option><option>Strong</option></select></div>
    <div><label>Case-control OR / CI lower</label><div style="display:flex;gap:4px"><input type="number" step="any" id="ce-or" placeholder="OR"><input type="number" step="any" id="ce-orlo" placeholder="CI low"></div></div></div>
   <div style="margin-top:6px">
    <label style="font-weight:400"><input type="checkbox" id="ce-pp4"> Phenotype highly specific for this gene (PP4)</label>
    <label style="font-weight:400"><input type="checkbox" id="ce-bs4"> Does not segregate with disease (BS4)</label>
    <label style="font-weight:400"><input type="checkbox" id="ce-bp2"> In trans with a pathogenic variant (dominant) or in cis with one (BP2)</label>
    <label style="font-weight:400"><input type="checkbox" id="ce-bp5"> Alternative molecular cause found (BP5)</label></div>
  </details>
  <label>Reviewer evidence <span class="hint">(advanced: ACMG codes directly) semicolon-separated, e.g. c.20A&gt;T=PS3,PP1_Strong; p.Glu7Val=-PM2; *=PS2</span></label><input type="text" id="m-evidence">
  <label>Patient sequence(s) * <span class="hint">DNA or protein; FASTA or raw; several sequences allowed; IUPAC (R, Y, W...) = heterozygous; end a complete protein with *</span></label>
  <textarea id="m-query" style="min-height:130px" placeholder=">Patient_exon1&#10;ATGGTGCATCTGACTCCTG..."></textarea>
  <input class="file" type="file" id="m-file" multiple accept=".fa,.fasta,.fna,.faa,.txt,.seq,.fas,.ab1,.abi">
  <div class="hint">Sequence files fill the box above; <b>Sanger trace files (.ab1)</b> are attached below - forward and/or reverse reads,
   heterozygous SNVs and indels are called from the peaks. The sequence box can stay empty when traces are attached.</div>
  <div id="m-traces" style="margin-top:6px"></div>
  <label>Notes <span class="hint">(optional, kept with the patient's inputs)</span></label><input type="text" id="m-notes">
  <div id="m-err" class="warn" style="margin-top:8px"></div>
 </div>
 <div class="mfoot"><button class="b secondary" id="m-cancel">Cancel</button>
  <button class="b secondary" id="m-save-next">Save &amp; add another</button><button class="b primary" id="m-save">Save patient</button></div>
</div></div>

<div class="modal-bg" id="dup-modal"><div class="modal" style="width:min(640px,94vw)">
 <div class="mhead"><span>Some patients are already in the list</span></div>
 <div class="mbody">
  <p id="dup-text" style="font-size:14px"></p>
  <div id="dup-list" class="mono" style="max-height:150px;overflow:auto;background:#f6f8fb;border-radius:8px;padding:8px 12px;font-size:12.5px"></div>
  <div id="dup-fam" style="margin-top:12px;display:none"><label style="font-weight:400"><input type="checkbox" id="dup-fam-rename" checked>
   <b>Family ID(s) <span id="dup-fam-ids"></span> are also already used</b> - these are different families: give the new ones a new Family ID
   <span class="hint">(untick only if they really are the same family)</span></label></div>
  <p class="hint" style="margin-top:12px">Patient IDs must be unique in a run. <b>Rename</b> adds the family ID (or a number) in front of each new
   duplicate ID, e.g. <code>P1_Father</code> &rarr; <code>FAM-HBB-01_P1_Father</code>.</p>
 </div>
 <div class="mfoot"><button class="b secondary" data-dup="cancel">Cancel import</button><button class="b secondary" data-dup="skip">Skip duplicates</button>
  <button class="b secondary" data-dup="replace">Replace existing</button><button class="b primary" data-dup="rename">Rename new patients</button></div>
</div></div>

<div class="view" id="v-results">
 <div class="toolbar"><div class="top"><div class="title" id="r-title">Results</div>
  <label style="margin:0;font-size:13px">Patient&nbsp;<select id="r-patient" style="width:auto;min-width:200px"></select></label>
  <button class="b secondary small" id="open-folder">Open results folder</button>
  <button class="b secondary small" id="open-browser">Print / save PDF</button></div>
  <div class="chips" id="r-chips"></div>
  <div class="tabs"><button data-p="cohort">Cohort overview</button><button data-p="maps" class="on">Variant maps</button><button data-p="report">Full report</button>
   <button data-p="aln">Alignments</button><button data-p="summary">Clinical summary</button><button data-p="files">Files</button></div></div>
 <div class="pane" id="p-cohort"></div>
 <div class="pane on" id="p-maps"></div>
 <div class="pane" id="p-report"><iframe id="report-frame"></iframe></div>
 <div class="pane" id="p-aln"></div>
 <div class="pane" id="p-summary"></div>
 <div class="pane" id="p-files"></div>
</div>

<div class="view" id="v-history"><div class="wrap"><div class="card"><h2>Previous analyses</h2>
<p class="hint" id="h-folder"></p><table class="hist"><thead><tr><th>Date</th><th>Patients</th><th>Variants</th><th>Highest class</th><th>Results</th></tr></thead>
<tbody id="h-body"></tbody></table></div></div></div>

<div class="view" id="v-help"><div class="wrap"><div class="card help"><h2>How to use</h2><ol>
<li><b>Reference genes</b>: add one reference per gene - FASTA file(s), an NM_ accession, or pasted sequence. Several genes can be analysed in one run; each patient is matched to the right gene automatically (or choose the gene in the patient pop-up). Set the gene symbol, inheritance (e.g. AD/AR for SOD1 p.Asp91Ala) and CDS per reference. <i>Save to library</i> keeps references for future sessions. Try <i>Load multi-gene example</i>.</li>
<li><b>Automatic annotation</b> (on by default) finds the gene's known pathogenic variants (ClinVar), transcript structure (RefSeq), domains and diseases (UniProt), constraint (gnomAD) and, for each detected variant, gnomAD frequency, REVEL, SpliceAI and AlphaMissense (Ensembl VEP). Use <i>Preview online annotation</i> to check what was found.</li>
<li>Add your own known pathogenic variants in the manual box (or load a .txt/.csv/.tsv file). An annotation JSON file can still override any value.</li>
<li><b>Clinical evidence</b>: in the patient window, open <i>Clinical evidence (guided)</i> and enter what you observed (meioses, cases in trans, functional assay, de novo, case data...). CLAVAGEN converts it to ACMG criteria with the ClinGen rules. A <b>gene rule profile</b> (section 3) applies an expert panel's gene-specific settings.</li>
<li><b>Patients</b>: click <i>+ Add patient</i> - a pop-up window asks for the patient's ID, sex, age, indication, phenotypes, reviewer evidence and sequence(s). Add as many patients as needed; use <i>Edit</i>, <i>Copy</i> or <i>Remove</i> in the list. <i>Import patient list</i> reads a CSV/TSV with columns <code>id, sex, age, indication, phenotypes, evidence, query</code> (query = sequence).</li>
<li>Choose where results are saved with <i>Choose folder...</i> (remembered), then <b>Run analysis</b>.</li>
<li><b>Results</b>: pick a patient in the <i>Patient</i> list to see their variant maps, full report, alignments, clinical summary and files. <i>Cohort overview</i> shows all patients together.</li>
<li>Each run is saved in one dated folder with a sub-folder per patient, plus <code>cohort_summary.html</code> / <code>.tsv</code>. <b>History</b> re-opens any previous run.</li></ol>
<p class="hint">Heterozygous SNVs are recognised from IUPAC codes (R, Y, S, W, K, M). Results are decision support only and must be reviewed by a qualified clinical scientist / geneticist.</p></div></div></div>

<script>
const $=id=>document.getElementById(id);
const COLORS={"Pathogenic":"#b2182b","Likely pathogenic":"#ef6548","Uncertain significance":"#e3a008","Likely benign":"#6baed6","Benign":"#2166ac"};
const esc=s=>String(s??"").replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const badge=c=>'<span class="badge" style="background:'+(COLORS[c]||"#999")+'">'+esc(c)+'</span>';
const stars=n=>n?"\u2605".repeat(n):"0\u2605";
const basisTag=v=>(v.basis&&!v.basis.startsWith("ACMG"))?'<div class="hint">'+(v.basis.startsWith("Known-variant")?'from your known-variant list':'from ClinVar '+stars(v.stars)+((v.stars||0)<2?' (limited review)':''))+' &middot; ACMG auto: '+esc(v.acmg)+'</div>':
 (v.clinvar?'<div class="hint">ClinVar: '+esc(v.clinvar)+' '+stars(v.stars)+'</div>':'');
const post=(u,b)=>fetch(u,{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(b)}).then(r=>r.json());
let current=null, cur=null;
function show(v){document.querySelectorAll('.nav button').forEach(b=>b.classList.toggle('on',b.dataset.v===v));
 document.querySelectorAll('.view').forEach(x=>x.classList.toggle('on',x.id==='v-'+v)); if(v==='history')loadHistory();}
document.querySelectorAll('.nav button').forEach(b=>b.onclick=()=>{if(!b.disabled)show(b.dataset.v)});
function pane(p){document.querySelectorAll('.tabs button').forEach(b=>b.classList.toggle('on',b.dataset.p===p));
 document.querySelectorAll('.pane').forEach(x=>x.classList.toggle('on',x.id==='p-'+p));
 if(p==='report'&&cur&&cur.report&&!$("report-frame").src.endsWith(cur.report))$("report-frame").src=cur.report;}
document.querySelectorAll('.tabs button').forEach(b=>b.onclick=()=>pane(b.dataset.p));
document.querySelectorAll('input[type=file][data-target]').forEach(f=>f.onchange=()=>{const r=new FileReader();r.onload=()=>{$(f.dataset.target).value=r.result};if(f.files[0])r.readAsText(f.files[0]);});

/* ---------------- patients ---------------- */
let patients=[];try{patients=JSON.parse(sessionStorage.getItem("acmg_patients")||"[]")}catch(e){}
let editIdx=-1;
function savePatients(){try{sessionStorage.setItem("acmg_patients",JSON.stringify(patients))}catch(e){}}
function seqInfo(t,p){if(p&&p.traces&&p.traces.length)return p.traces.length+" trace(s)"+(t?" + "+seqInfo(t):"");if(!t)return "-";const recs=t.trim().startsWith(">")?t.split(/\n(?=>)/):[t];let n=0,len=0;
 for(const r of recs){const s=r.replace(/^>.*$/m,"").replace(/[^A-Za-z*]/g,"");if(s){n++;len+=s.length}}return n+" seq, "+len+" residues";}
function renderPatients(){
 $("pcount").textContent=patients.length?("- "+patients.length+" patient"+(patients.length>1?"s":"")):"";
 $("patient-rows").innerHTML=patients.length?patients.map((p,i)=>'<tr><td>'+(i+1)+'</td><td><b>'+esc(p.id)+'</b></td><td>'+(p.family?esc(p.family)+'<br><span class="hint">'+esc(p.relationship||'')+(p.affected&&p.affected!=="unknown"?', '+esc(p.affected):'')+'</span>':'')+'</td><td>'+esc(p.sex)+'</td><td>'+esc(p.age)+
  '</td><td>'+esc(p.indication)+(p.phenotypes?'<br><span class="hint">'+esc(p.phenotypes)+'</span>':'')+'</td><td class="hint">'+esc(seqInfo(p.query,p))+(p.reference?'<br>&rarr; '+esc((refs.find(r=>r.id===p.reference)||{gene:p.reference}).gene||p.reference):'')+
  '</td><td class="hint">'+esc(p.evidence)+(p.clinical_evidence?(p.evidence?'<br>':'')+esc(ceSummary(p.clinical_evidence)):'')+'</td><td style="white-space:nowrap"><button class="lnk" data-e="'+i+'">Edit</button> <button class="lnk" data-c="'+i+
  '">Copy</button> <button class="lnk red" data-d="'+i+'">Remove</button></td></tr>').join(''):
  '<tr><td colspan=9 class="empty" style="padding:26px">No patients yet - click <b>+ Add patient</b>.</td></tr>';
 document.querySelectorAll('[data-e]').forEach(b=>b.onclick=()=>openModal(+b.dataset.e));
 document.querySelectorAll('[data-c]').forEach(b=>b.onclick=()=>{const p=Object.assign({},patients[+b.dataset.c]);
   p.id=p.id+"_copy";patients.splice(+b.dataset.c+1,0,p);savePatients();renderPatients();openModal(+b.dataset.c+1);});
 document.querySelectorAll('[data-d]').forEach(b=>b.onclick=()=>{const i=+b.dataset.d;
   if(confirmRemove(patients[i].id)){patients.splice(i,1);savePatients();renderPatients();}});
}
function confirmRemove(id){return window.confirm?window.confirm("Remove patient "+id+"?"):true;}
const MF=["id","sex","age","family","relationship","affected","reference","indication","phenotypes","evidence","query","notes"];
function openModal(i){editIdx=i;const p=i>=0?patients[i]:{};
 $("m-title").textContent=i>=0?"Edit patient":"Add patient";$("m-err").textContent="";
 MF.forEach(k=>$("m-"+k).value=p[k]||(k==="affected"?"unknown":""));$("m-parentage").checked=!!p.parentage_confirmed;$("m-file").value="";
 mTraces=(p.traces||[]).slice();renderTraces();fillCE(p.clinical_evidence);
 $("m-reference").innerHTML='<option value="">Auto-detect (best-matching reference)</option>'+refs.map(r=>'<option value="'+esc(r.id)+'">'+esc((r.gene||'?')+' - '+r.id)+'</option>').join('');
 $("m-reference").value=p.reference||"";
 $("fam-list").innerHTML=[...new Set(patients.map(x=>x.family).filter(Boolean))].map(f=>'<option value="'+esc(f)+'">').join('');
 $("modal").classList.add("on");setTimeout(()=>$("m-id").focus(),50);}
function closeModal(){$("modal").classList.remove("on");}
const CEN=["meioses","families","healthy","tp-c","tp-u","tlp-c","tlp-u","hom","oddspath","dn-c","dn-a","ap","or","orlo"];
function num(id){const v=$("ce-"+id).value;return v===""?null:+v;}
function readCE(){const e={};
 if(num("meioses")!==null)e.segregation={informative_meioses:num("meioses"),families:num("families")||undefined};
 const tr=[];[["tp-c","Pathogenic","confirmed"],["tp-u","Pathogenic","unknown"],["tlp-c","Likely pathogenic","confirmed"],["tlp-u","Likely pathogenic","unknown"]]
  .forEach(([id,pt,ph])=>{for(let k=0;k<(num(id)||0);k++)tr.push({partner:pt,phase:ph});});
 if(tr.length)e.in_trans=tr; if(num("hom"))e.homozygous_probands=num("hom");
 if($("ce-func").value){e.functional={result:$("ce-func").value};if(num("oddspath")!==null)e.functional.oddspath=num("oddspath");else if($("ce-fstr").value)e.functional.strength=$("ce-fstr").value;}
 const dn=[];for(let k=0;k<(num("dn-c")||0);k++)dn.push({confirmed:true,phenotype:$("ce-dn-ph").value});
 for(let k=0;k<(num("dn-a")||0);k++)dn.push({confirmed:false,phenotype:$("ce-dn-ph").value});if(dn.length)e.de_novo=dn;
 if(num("ap"))e.affected_probands={count:num("ap"),strength:$("ce-apstr").value};
 if(num("or")!==null)e.case_control={odds_ratio:num("or"),ci_lower:num("orlo")||0};
 if(num("healthy"))e.healthy_adults={count:num("healthy")};
 if($("ce-pp4").checked)e.phenotype_specific=true;if($("ce-bs4").checked)e.non_segregation=true;
 if($("ce-bp2").checked)e.benign_phase=true;if($("ce-bp5").checked)e.alternative_cause=true;
 if(!Object.keys(e).length)return null;const o={};o[$("ce-key").value.trim()||"*"]=e;return o;}
function fillCE(ce){CEN.forEach(id=>$("ce-"+id).value="");["func","fstr"].forEach(id=>$("ce-"+id).value="");$("ce-apstr").value="Supporting";
 ["pp4","bs4","bp2","bp5"].forEach(id=>$("ce-"+id).checked=false);$("ce-key").value="";$("ce-dn-ph").value="specific";
 if(!ce||!Object.keys(ce).length){$("m-ce-box").open=false;return;}
 const k=Object.keys(ce)[0],e=ce[k]||{};$("ce-key").value=k==="*"?"":k;$("m-ce-box").open=true;
 if(e.segregation){$("ce-meioses").value=e.segregation.informative_meioses??"";$("ce-families").value=e.segregation.families??"";}
 const cnt=(pt,ph)=>(e.in_trans||[]).filter(t=>t.partner===pt&&t.phase===ph).length||"";
 $("ce-tp-c").value=cnt("Pathogenic","confirmed");$("ce-tp-u").value=cnt("Pathogenic","unknown");
 $("ce-tlp-c").value=cnt("Likely pathogenic","confirmed");$("ce-tlp-u").value=cnt("Likely pathogenic","unknown");
 $("ce-hom").value=e.homozygous_probands||"";
 if(e.functional){$("ce-func").value=e.functional.result||"";$("ce-oddspath").value=e.functional.oddspath??"";$("ce-fstr").value=e.functional.strength||"";}
 if(e.de_novo){$("ce-dn-c").value=e.de_novo.filter(d=>d.confirmed).length||"";$("ce-dn-a").value=e.de_novo.filter(d=>!d.confirmed).length||"";$("ce-dn-ph").value=(e.de_novo[0]||{}).phenotype||"specific";}
 if(e.affected_probands){$("ce-ap").value=e.affected_probands.count;$("ce-apstr").value=e.affected_probands.strength||"Supporting";}
 if(e.case_control){$("ce-or").value=e.case_control.odds_ratio;$("ce-orlo").value=e.case_control.ci_lower;}
 if(e.healthy_adults)$("ce-healthy").value=e.healthy_adults.count;
 $("ce-pp4").checked=!!e.phenotype_specific;$("ce-bs4").checked=!!e.non_segregation;$("ce-bp2").checked=!!e.benign_phase;$("ce-bp5").checked=!!e.alternative_cause;}
function ceSummary(ce){if(!ce)return "";const e=Object.values(ce)[0]||{},k=Object.keys(ce)[0],L=[];
 if(e.segregation)L.push(e.segregation.informative_meioses+" meioses");if(e.in_trans)L.push(e.in_trans.length+" in trans");
 if(e.homozygous_probands)L.push(e.homozygous_probands+" hom");if(e.functional)L.push("assay "+e.functional.result);
 if(e.de_novo)L.push(e.de_novo.length+" de novo");if(e.affected_probands)L.push(e.affected_probands.count+" probands");
 if(e.case_control)L.push("OR "+e.case_control.odds_ratio);if(e.healthy_adults)L.push(e.healthy_adults.count+" healthy");
 if(e.phenotype_specific)L.push("PP4");if(e.non_segregation)L.push("BS4");if(e.benign_phase)L.push("BP2");if(e.alternative_cause)L.push("BP5");
 return (k==="*"?"":k+": ")+L.join(", ");}
function modalSave(again){const p={};MF.forEach(k=>p[k]=$("m-"+k).value.trim());p.parentage_confirmed=$("m-parentage").checked;
 if(p.family&&!p.relationship){$("m-err").textContent="Please choose the relationship to the proband for family members.";return;}
 if(!p.id){$("m-err").textContent="Patient ID is required.";return;}
 if(patients.some((x,j)=>x.id===p.id&&j!==editIdx)){$("m-err").textContent="Patient ID '"+p.id+"' is already in the list.";return;}
 p.traces=mTraces.slice();p.clinical_evidence=readCE();
 if(!p.query&&!p.traces.length){$("m-err").textContent="Please paste or load the patient's sequence, or attach a Sanger trace (.ab1).";return;}
 if(editIdx>=0)patients[editIdx]=p;else patients.push(p);savePatients();renderPatients();
 if(again)openModal(-1);else closeModal();}
$("add-patient").onclick=()=>openModal(-1);
$("m-close").onclick=closeModal;$("m-cancel").onclick=closeModal;
$("m-save").onclick=()=>modalSave(false);$("m-save-next").onclick=()=>modalSave(true);
$("modal").onclick=e=>{if(e.target.id==="modal")closeModal();};
document.addEventListener("keydown",e=>{if(e.key==="Escape"&&$("modal").classList.contains("on"))closeModal();});
let mTraces=[];
function renderTraces(){$("m-traces").innerHTML=mTraces.map((t,i)=>'<span class="chip">&#127897; '+esc(t.name)+' ('+Math.round(t.data.length*0.75/1024)+' KB) <button class="lnk red" data-tr="'+i+'">remove</button></span>').join(' ');
 document.querySelectorAll('[data-tr]').forEach(b=>b.onclick=()=>{mTraces.splice(+b.dataset.tr,1);renderTraces();});}
$("m-file").onchange=()=>{for(const f of $("m-file").files){const r=new FileReader();
  if(/\.(ab1|abi)$/i.test(f.name)){r.onload=()=>{mTraces.push({name:f.name,data:r.result});renderTraces();
    if(!$("m-id").value.trim())$("m-id").value=f.name.replace(/[._-]?(F|R|fwd|rev|forward|reverse)?\.(ab1|abi)$/i,"");};r.readAsDataURL(f);}
  else{r.onload=()=>{$("m-query").value=r.result;if(!$("m-id").value.trim())$("m-id").value=f.name.replace(/\.[^.]+$/,"");};r.readAsText(f);}}};
$("clear-patients").onclick=()=>{if(patients.length&&confirmRemove("ALL ("+patients.length+")")){patients=[];savePatients();renderPatients();}};
$("import-btn").onclick=()=>$("import-file").click();
function askDuplicates(name,dups,famDups){return new Promise(res=>{
 $("dup-text").innerHTML='<b>'+dups.length+'</b> patient ID(s) in <b>'+esc(name)+'</b> already exist in the patient list:';
 $("dup-list").innerHTML=dups.map(esc).join('<br>');
 $("dup-fam").style.display=famDups.length?'':'none';$("dup-fam-ids").textContent=famDups.join(', ');$("dup-fam-rename").checked=true;
 $("dup-modal").classList.add("on");
 document.querySelectorAll('[data-dup]').forEach(b=>b.onclick=()=>{$("dup-modal").classList.remove("on");res({action:b.dataset.dup,renameFam:$("dup-fam-rename").checked});});});}
function uniqueId(base){let id=base,k=2;while(patients.some(x=>x.id===id)){id=base+"_"+k;k++;}return id;}
async function importPatients(text,name){
 const d=await post("/parsepatients",{text:text});
 if(!d.ok)return {error:d.error};
 let incoming=d.patients;
 const existingIds=new Set(patients.map(x=>x.id)), existingFams=new Set(patients.map(x=>x.family).filter(Boolean));
 const dups=incoming.filter(p=>existingIds.has(p.id)).map(p=>p.id);
 const famDups=[...new Set(incoming.map(p=>p.family).filter(f=>f&&existingFams.has(f)))];
 let renamed=[],replaced=0,skipped=0,famRenamed=[];
 if(dups.length||famDups.length){
  const ch=await askDuplicates(name,dups,famDups);
  if(ch.action==="cancel")return {cancelled:true};
  if(famDups.length&&ch.renameFam&&ch.action!=="replace"){
   const map={};for(const f of famDups){let nf=f+"_2",k=3;while(existingFams.has(nf)){nf=f+"_"+k;k++;}map[f]=nf;famRenamed.push(f+" → "+nf);}
   incoming=incoming.map(p=>map[p.family]?Object.assign({},p,{family:map[p.family]}):p);}
  if(ch.action==="skip"){skipped=incoming.filter(p=>existingIds.has(p.id)).length;incoming=incoming.filter(p=>!existingIds.has(p.id));}
  else if(ch.action==="replace"){for(const p of incoming){const i=patients.findIndex(x=>x.id===p.id);if(i>=0){patients.splice(i,1);replaced++;}}}
  else{incoming=incoming.map(p=>{if(!existingIds.has(p.id))return p;
    const nid=uniqueId((p.family?p.family+"_":"")+p.id);renamed.push(p.id+" → "+nid);return Object.assign({},p,{id:nid});});}
 }
 for(const p of incoming){if(patients.some(x=>x.id===p.id)){const nid=uniqueId(p.id);renamed.push(p.id+" → "+nid);p.id=nid;}patients.push(p);}
 savePatients();renderPatients();
 return {added:incoming.length,renamed,replaced,skipped,famRenamed};}
$("import-file").onchange=async()=>{const files=[...$("import-file").files];if(!files.length)return;const msgs=[];
 for(const f of files){const r=await importPatients(await f.text(),f.name);
  if(r.error){msgs.push(f.name+': import failed - '+r.error);continue;}
  if(r.cancelled){msgs.push(f.name+': import cancelled');continue;}
  let m=f.name+': imported '+r.added+' patient(s)';
  if(r.renamed.length)m+='; renamed '+r.renamed.length+' duplicate ID(s) ('+r.renamed.slice(0,3).join(', ')+(r.renamed.length>3?', ...':'')+')';
  if(r.replaced)m+='; replaced '+r.replaced+' existing';
  if(r.skipped)m+='; skipped '+r.skipped+' duplicate(s)';
  if(r.famRenamed.length)m+='; new Family ID(s): '+r.famRenamed.join(', ');
  msgs.push(m+'.');}
 $("status").innerHTML=msgs.map(esc).join('<br>')+'<br>'+patients.length+' patient(s) in the list.';$("import-file").value="";};
renderPatients();

/* ---------------- references ---------------- */
let refs=[];try{refs=JSON.parse(sessionStorage.getItem("acmg_refs")||"[]")}catch(e){}
function saveRefs(){try{sessionStorage.setItem("acmg_refs",JSON.stringify(refs))}catch(e){}}
const INH=[["","from databases"],["AD","AD - dominant"],["AR","AR - recessive"],["AD/AR","AD/AR - both"],["XLR","X-linked recessive"],["XLD","X-linked dominant"]];
function renderRefs(){
 $("rcount").textContent=refs.length?("- "+refs.length+" reference"+(refs.length>1?"s":"")):"";
 $("ref-rows").innerHTML=refs.length?refs.map((r,i)=>'<tr><td><b>'+esc(r.id)+'</b><br><span class="hint">'+esc((r.description||'').slice(0,60))+'</span></td>'+
  '<td><input type="text" data-rf="gene" data-ri="'+i+'" value="'+esc(r.gene||'')+'" placeholder="symbol" style="width:90px"></td>'+
  '<td class="hint">'+esc(r.accession?'accession':(r.type+', '+r.length+(r.type==='protein'?' aa':' nt')))+'</td>'+
  '<td><select data-rf="inheritance" data-ri="'+i+'" style="width:auto">'+INH.map(o=>'<option value="'+o[0]+'"'+((r.inheritance||'')===o[0]?' selected':'')+'>'+o[1]+'</option>').join('')+'</select></td>'+
  '<td><input type="text" data-rf="cds" data-ri="'+i+'" value="'+esc(r.cds||'')+'" placeholder="auto" style="width:90px"></td>'+
  '<td style="white-space:nowrap"><button class="lnk" data-rp="'+i+'">Preview</button> <button class="lnk red" data-rd="'+i+'">Remove</button></td></tr>').join(''):
  '<tr><td colspan=6 class="empty" style="padding:22px">No reference yet - add a FASTA file, an accession (e.g. NM_000518.5) or paste a sequence.</td></tr>';
 document.querySelectorAll('[data-rf]').forEach(el=>el.onchange=()=>{refs[+el.dataset.ri][el.dataset.rf]=el.value.trim();saveRefs();});
 document.querySelectorAll('[data-rd]').forEach(b=>b.onclick=()=>{refs.splice(+b.dataset.rd,1);saveRefs();renderRefs();renderPatients();});
 document.querySelectorAll('[data-rp]').forEach(b=>b.onclick=()=>previewRefs([refs[+b.dataset.rp]]));}
async function addRefs(text,name){const r=await post("/parsereferences",{text:text,name:name});
 if(!r.ok){$("ref-msg").innerHTML='<span class="warn">'+esc(r.error)+'</span>';return;}
 let added=0,skipped=[];for(const x of r.references){if(refs.some(y=>y.id===x.id)){skipped.push(x.id);continue;}refs.push(x);added++;}
 saveRefs();renderRefs();$("ref-msg").textContent="Added "+added+" reference(s)"+(skipped.length?"; already in the list: "+skipped.join(", "):"")+".";return added;}
$("ref-add-file").onclick=()=>$("ref-file").click();
$("ref-file").onchange=async()=>{let n=0;const files=[...$("ref-file").files];for(const f of files){n+=(await addRefs(await f.text(),f.name.replace(/\.[^.]+$/,"")))||0;}
 $("ref-msg").textContent="Added "+n+" reference(s) from "+files.length+" file(s). Check the gene and set inheritance if needed.";$("ref-file").value="";};
$("ref-add-acc").onclick=()=>{const v=$("ref-acc").value.trim();if(v){addRefs(v.split(/[\s,;]+/).join("\n"));$("ref-acc").value="";}};
$("ref-paste-toggle").onclick=()=>{$("ref-paste").style.display=$("ref-paste").style.display==="none"?"":"none";};
$("ref-add-paste").onclick=()=>{if($("reference").value.trim()){addRefs($("reference").value);$("reference").value="";}};
let library=[];
async function loadLibrary(){const r=await post("/library",{});library=r.library||[];
 $("ref-lib").innerHTML='<option value="">Add from library'+(library.length?' ('+library.length+')':' (empty)')+'...</option>'+library.map((x,i)=>'<option value="'+i+'">'+esc((x.gene||'?')+' - '+x.id)+'</option>').join('');}
$("ref-lib").onchange=()=>{const i=$("ref-lib").value;if(i==="")return;const x=library[+i];
 if(!refs.some(y=>y.id===x.id)){refs.push(Object.assign({},x));saveRefs();renderRefs();$("ref-msg").textContent="Added "+x.id+" from the library.";}$("ref-lib").value="";};
$("ref-save-lib").onclick=async()=>{if(!refs.length){$("ref-msg").textContent="Nothing to save.";return;}
 await post("/library",{action:"add",references:refs});await loadLibrary();$("ref-msg").textContent="Saved "+refs.length+" reference(s) to your library (available in future sessions).";};
renderRefs();loadLibrary();

/* ---------------- settings / annotation ---------------- */
const fields=["annotation","strict","known","rules"];
$("clear").onclick=()=>{fields.forEach(k=>$(k).value=k==="strict"?"0":"");patients=[];savePatients();renderPatients();refs=[];saveRefs();renderRefs();};
async function loadSettings(){const st=await (await fetch("/settings")).json();$("results_dir").value=st.results_dir;}
loadSettings();
$("choose").onclick=async()=>{$("status").textContent="Choose a folder in the dialog...";
 const r=await post("/choosefolder",{});
 if(r.ok&&r.path){$("results_dir").value=r.path;$("status").textContent="Results will be saved in "+r.path;}
 else $("status").textContent=r.error?("Folder dialog unavailable: "+r.error+" - type the path and click 'Use this folder'."):"No folder chosen.";};
$("setdir").onclick=async()=>{const r=await post("/settings",{results_dir:$("results_dir").value});
 $("status").textContent=r.ok?("Results will be saved in "+r.results_dir):("Error: "+r.error);if(r.ok)$("results_dir").value=r.results_dir;};
async function previewRefs(list){
 if(!list.length){$("pstatus").textContent="Add a reference first (section 1).";return;}
 $("pstatus").textContent="Searching ClinVar, RefSeq, UniProt and gnomAD for "+list.length+" reference(s)... (up to a minute the first time)";
 $("preview").disabled=true;const b=$("preview-box");b.style.display="";b.innerHTML="";
 try{for(const ref of list){
  const r=await post("/annotate",{gene:ref.gene||"",reference:ref.accession||('>'+ref.id+' '+(ref.description||'')),annotation:$("annotation").value});
  let h='<div style="margin-bottom:10px">';
  if(!r.ok){h+='<b>'+esc(ref.id)+'</b>: <span class="warn">'+esc(r.error)+'</span></div>';b.innerHTML+=h;continue;}
  if(!ref.gene&&r.gene){ref.gene=r.gene;saveRefs();renderRefs();}
  h+='<b>'+esc(r.gene)+'</b> '+esc(r.transcript||"")+' &middot; inheritance: '+esc(ref.inheritance||r.inheritance||"?")+(ref.inheritance?' (set in list)':'')+' &middot; LOF mechanism: '+esc(r.lof)+'<br>';
  h+=r.diseases.map(d=>'&bull; '+esc(d)).join('<br>')+(r.diseases.length?'<br>':'');
  h+='<ul style="margin:6px 0">'+r.sources.map(x=>'<li>'+esc(x)+'</li>').join('')+r.errors.map(x=>'<li class="warn">Not available: '+esc(x)+'</li>').join('')+'</ul>';
  h+='<b>'+r.n_known+'</b> known P/LP variants found'+(r.n_known?' (first 25):':'')+'<div style="max-height:130px;overflow:auto;font-family:Menlo,Consolas,monospace;font-size:12px">'+
   r.known.map(k=>esc(k.hgvs_c+'  '+(k.hgvs_p||'')+'  '+k.classification+'  '+k.source)).join('<br>')+'</div></div>';
  b.innerHTML+=h;}
  $("pstatus").textContent="Done.";}catch(e){$("pstatus").textContent="Error: "+e;}finally{$("preview").disabled=false;}}
$("preview").onclick=()=>previewRefs(refs);
$("checkconn").onclick=async()=>{$("pstatus").textContent="Testing connections to NCBI, UniProt, gnomAD and Ensembl...";
 const r=await post("/checkconn",{});const b=$("preview-box");b.style.display="";
 b.innerHTML='<b>Connection check</b><ul style="margin:6px 0">'+r.results.map(x=>'<li class="'+(x[1]?'':'warn')+'">'+(x[1]?'&#10003; ':'&#10007; ')+
  esc(x[0])+': '+esc(x[2])+'</li>').join('')+'</ul>';$("pstatus").textContent="";};
$("example").onclick=async()=>{const d=await (await fetch("/example")).json();fields.forEach(k=>{if(d[k]!==undefined)$(k).value=d[k]});refs=d.references||[];saveRefs();renderRefs();
 patients=d.patients;savePatients();renderPatients();$("status").textContent="Example loaded ("+patients.length+" patients) - click Run analysis.";};
$("famexample").onclick=async()=>{const d=await (await fetch("/example?family=1")).json();fields.forEach(k=>{if(d[k]!==undefined)$(k).value=d[k]});refs=d.references||[];saveRefs();renderRefs();
 patients=d.patients;savePatients();renderPatients();$("status").textContent="Family example loaded (trio + sibling, and a carrier couple) - click Run analysis.";};
$("multiexample").onclick=async()=>{const d=await (await fetch("/example?multi=1")).json();fields.forEach(k=>{if(d[k]!==undefined)$(k).value=d[k]});
 refs=d.references||[];saveRefs();renderRefs();patients=d.patients;savePatients();renderPatients();
 $("status").textContent="Multi-gene example loaded: 2 references (HBB, SOD1 set to AD/AR) and 12 patients in 2 families - click Run analysis.";};
$("open-folder").onclick=()=>current&&post("/open",{what:"folder",run:current.name});
$("open-browser").onclick=()=>current&&cur&&post("/open",{what:"report",run:current.name,path:cur.base.split("/runs/"+current.name+"/")[1]+"report.html"});

/* ---------------- results ---------------- */
async function svgToPng(svgText){return new Promise((res,rej)=>{const img=new Image();
 img.onload=()=>{try{const s=2,c=document.createElement('canvas');c.width=img.naturalWidth*s;c.height=img.naturalHeight*s;
  const x=c.getContext('2d');x.fillStyle='#fff';x.fillRect(0,0,c.width,c.height);x.scale(s,s);x.drawImage(img,0,0);res(c.toDataURL('image/png'));}catch(e){rej(e)}};
 img.onerror=rej;img.src='data:image/svg+xml;base64,'+btoa(unescape(encodeURIComponent(svgText)));});}
const rel=u=>u.split("/runs/"+current.name+"/")[1];
function renderCohort(d){
 let h='<div class="wrap">'+(d.cohort?'<p><button class="b secondary small" id="open-cohort">Open cohort_summary.html</button></p>':'')+
  '<table><tr><th>Patient</th><th>Sex / age</th><th>Indication / phenotypes</th><th>Variants (ACMG class)</th><th>Result</th></tr>';
 d.patients.forEach((p,i)=>{const pt=p.patient||{};
  if(!p.ok){h+='<tr><td><b>'+esc(p.id)+'</b></td><td colspan=4 class="warn">Analysis failed: '+esc(p.error)+'</td></tr>';return;}
  const seen=new Set();const vs=p.samples.flatMap(s=>s.variants).filter(v=>{const k=(v.gene||'')+'|'+v.c+'|'+v.p;if(seen.has(k))return false;seen.add(k);return true;});
  h+='<tr><td><a href="#" data-pi="'+i+'"><b>'+esc(p.id)+'</b></a></td><td>'+esc([pt.sex,pt.age].filter(Boolean).join(', '))+'</td><td>'+esc(pt.indication||'')+
   '<br><span class="hint">'+esc((p.phenotypes||[]).join('; '))+'</span></td><td>'+(vs.length?vs.map(v=>'<b>'+esc(v.gene||'')+'</b> <span class="mono">'+esc(v.c+' '+v.p)+'</span> '+badge(v.cls)+basisTag(v)).join('<div style="height:6px"></div>'):'No variants detected')+
   '</td><td>'+esc([...new Set(p.samples.map(s=>s.status))].join('; '))+'</td></tr>';});
 h+='</table><p class="hint">Click a patient to open their results.</p>';
 for(const f of (d.families||[])){
  h+='<div class="family"><h3>Family '+esc(f.family)+' <span class="hint">(proband: '+esc(f.proband||'-')+')</span></h3>'+
   (f.warnings||[]).map(w=>'<div class="warn">'+esc(w)+'</div>').join('')+'<img src="'+f.pedigree_url+'" alt="pedigree">'+
   '<h4>Findings</h4><ul>'+((f.findings||[]).map(x=>'<li>'+esc(x)+'</li>').join('')||'<li>None.</li>')+'</ul>'+
   ((f.evidence||[]).length?'<h4>ACMG evidence added from family data</h4><ul>'+f.evidence.map(e=>'<li><b>'+esc(e.code+'_'+e.strength)+'</b>: '+esc(e.rationale)+'</li>').join('')+'</ul>':'')+
   (f.couples||[]).map(c=>'<h4>Couple '+esc(c.label)+'</h4><ul>'+c.risks.map(r=>'<li class="'+(r.level==="high"?'risk':'')+'"><b>'+esc(r.risk)+'</b> - '+esc(r.text)+'</li>').join('')+'</ul>').join('');
  const ids=f.members.map(m=>m.id),sym={"het":"&#9680; het","hom":"&#9679; hom","absent":"&#9675; absent","carrier":"&#9679; carrier","n/c":"n.t."};
  h+='<h4>Segregation table</h4><table><tr><th>Variant</th>'+f.members.map(m=>'<th>'+esc(m.id)+'<br><span class="hint">'+esc(m.relationship)+', '+esc(m.affected)+'</span></th>').join('')+'</tr>'+
   f.variants.map(v=>'<tr><td><b>'+esc(v.gene)+'</b> <span class="mono">'+esc(v.hgvs_c+' '+v.hgvs_p)+'</span><br>'+badge(v.classification)+'</td>'+ids.map(i=>'<td>'+(sym[v.genotypes[i]]||esc(v.genotypes[i]||''))+'</td>').join('')+'</tr>').join('')+'</table></div>';}
 $("p-cohort").innerHTML=h+'</div>';
 document.querySelectorAll('[data-pi]').forEach(a=>a.onclick=e=>{e.preventDefault();$("r-patient").value=a.dataset.pi;renderPatient(+a.dataset.pi);pane('maps');});
 const oc=$("open-cohort");if(oc)oc.onclick=()=>post("/open",{what:"file",run:d.name,name:"cohort_summary.html"});}
async function renderPatient(i){
 const p=current.patients[i];cur=p;$("report-frame").src="about:blank";
 if(!p.ok){$("p-maps").innerHTML='<div class="wrap warn">Analysis failed for '+esc(p.id)+': '+esc(p.error)+'</div>';return;}
 $("r-chips").innerHTML=p.samples.map(s=>'<span class="chip"><b>'+esc(p.id)+' / '+esc(s.query)+'</b>: '+esc(s.status)+' '+(s.highest&&s.highest!=="None"?badge(s.highest):'')+'</span>').join('');
 let maps='',alns='';
 for(const s of p.samples){
  maps+='<div class="mapbox"><h3>'+esc(p.id)+': '+esc(s.query)+' vs '+esc(s.reference)+' &mdash; '+esc(s.status)+'</h3>'+
   s.warnings.map(w=>'<div class="warn">! '+esc(w)+'</div>').join('')+
   '<img src="'+s.map+'" alt="variant map"><table style="margin-top:8px"><tr><th>ID</th><th>c.</th><th>p.</th><th>Type</th><th>Zygosity</th><th>Class</th><th>ACMG criteria</th></tr>'+
   (s.variants.length?s.variants.map(v=>'<tr><td>'+esc(v.id)+'</td><td>'+esc(v.c)+'</td><td>'+esc(v.p)+'</td><td>'+esc(v.type)+'</td><td>'+esc(v.zyg)+'</td><td>'+badge(v.cls)+basisTag(v)+'</td><td>'+esc(v.crit)+'</td></tr>').join(''):'<tr><td colspan=7>No variants detected.</td></tr>')+'</table></div>';
  const t=await (await fetch(s.aln)).text();alns+='<pre class="aln">'+esc(t)+'</pre>';}
 $("p-maps").innerHTML=maps;$("p-aln").innerHTML=alns;
 const md=await (await fetch(p.base+'clinical_summary.md')).text();$("p-summary").innerHTML='<pre class="md">'+esc(md)+'</pre>';
 const f=await (await fetch('/files?run='+encodeURIComponent(current.name))).json();
 const pre=rel(p.base);const mine=f.files.filter(x=>x.name.startsWith(pre)||!x.name.includes("/"));
 $("p-files").innerHTML='<div class="wrap"><p>Folder: <code>'+esc(f.folder)+'</code></p><table><tr><th>File</th><th>Size</th></tr>'+
  mine.map(x=>'<tr><td><a href="#" data-f="'+esc(x.name)+'">'+esc(x.name)+'</a></td><td>'+x.size+'</td></tr>').join('')+'</table><p class="hint">Click a file to open it with its default program. Other patients\' files are in their own sub-folders.</p></div>';
 document.querySelectorAll('#p-files a[data-f]').forEach(a=>a.onclick=e=>{e.preventDefault();post("/open",{what:"file",run:current.name,name:a.dataset.f});});
 if($("p-report").classList.contains("on"))$("report-frame").src=p.report;
}
async function openRun(d,makePng){
 current=d;$("nav-results").disabled=false;show('results');
 $("r-title").innerHTML=esc(d.title);
 $("r-patient").innerHTML=d.patients.map((p,i)=>'<option value="'+i+'">'+esc(p.id)+(p.highest?' - '+esc(p.highest):(p.ok?'':' - failed'))+'</option>').join('');
 $("r-patient").onchange=()=>renderPatient(+$("r-patient").value);
 document.querySelector('.tabs button[data-p=cohort]').style.display=(d.patients.length>1||(d.families||[]).length)?'':'none';
 renderCohort(d);await renderPatient(0);pane(d.patients.length>1?'cohort':'maps');
 if(makePng){for(const p of d.patients){for(const s of p.samples){try{const svg=await (await fetch(s.map)).text();const png=await svgToPng(svg);
   await post("/savepng",{run:d.name,file:rel(s.map).replace(/\.svg$/,'.png'),data:png});}catch(e){console.log(e)}}}}}
$("run").onclick=async()=>{const body={};fields.forEach(k=>body[k]=$(k).value);body.online=$("online").checked?"1":"0";
 body.use_clinvar=$("clinvar_min").value==="off"?"0":"1";body.clinvar_min_stars=$("clinvar_min").value==="off"?"0":$("clinvar_min").value;
 body.results_dir=$("results_dir").value;body.patients=patients;body.references=refs;
 if(!refs.length){$("status").textContent="Please add at least one reference (section 1).";return;}
 if(!patients.length){$("status").textContent="Please add at least one patient (section 4).";openModal(-1);return;}
 $("status").textContent="Analysing "+patients.length+" patient(s) against "+refs.length+" reference(s)"+(body.online==="1"?" - searching online databases (up to a minute the first time)...":"...");$("run").disabled=true;
 try{const d=await post("/run",body);if(!d.ok){$("status").innerHTML='<span class="warn">Error: '+esc(d.error)+'</span>';return;}
  $("status").textContent="Done.";await openRun(d,true);}catch(e){$("status").textContent="Error: "+e;}finally{$("run").disabled=false;}};
async function loadHistory(){const h=await (await fetch("/history")).json();$("h-folder").textContent="Saved in: "+h.folder;
 $("h-body").innerHTML=h.runs.length?h.runs.map((r,i)=>'<tr data-i="'+i+'"><td>'+esc(r.date)+'</td><td>'+esc(r.patient_ids)+'</td><td>'+r.n_variants+'</td><td>'+(r.highest?badge(r.highest):'-')+
  '</td><td class="hint">'+esc(r.patients.map(p=>p.id+': '+(p.ok?[...new Set(p.samples.map(s=>s.status))].join('/'):'failed')).join(' | '))+'</td></tr>').join(''):'<tr><td colspan=5 class="empty">No analyses yet.</td></tr>';
 document.querySelectorAll('#h-body tr[data-i]').forEach(tr=>tr.onclick=()=>openRun(h.runs[+tr.dataset.i],false));}
</script></body></html>"""

from .paths import app_path  # noqa: E402
SETTINGS = app_path("settings.json")


LIBRARY = app_path("reference_library.json")


def load_library():
    try:
        with open(LIBRARY) as fh:
            return json.load(fh)
    except Exception:
        return []


def save_library(items):
    os.makedirs(os.path.dirname(LIBRARY), exist_ok=True)
    with open(LIBRARY, "w") as fh:
        json.dump(items, fh, indent=2)


def default_results_dir():
    docs = os.path.join(os.path.expanduser("~"), "Documents")
    return os.path.join(docs if os.path.isdir(docs) else os.path.expanduser("~"), "CLAVAGEN Results")


def load_settings():
    try:
        with open(SETTINGS) as fh:
            return json.load(fh)
    except Exception:
        return {}


def save_settings(st):
    os.makedirs(os.path.dirname(SETTINGS), exist_ok=True)
    with open(SETTINGS, "w") as fh:
        json.dump(st, fh, indent=2)


def set_results_dir(path):
    path = os.path.abspath(os.path.expanduser((path or "").strip()))
    if not path or path == os.path.abspath(""):
        raise ValueError("Please give a folder path.")
    os.makedirs(path, exist_ok=True)
    if not os.access(path, os.W_OK):
        raise ValueError(f"Cannot write to {path}")
    st = load_settings()
    st["results_dir"] = path
    save_settings(st)
    return path


def choose_folder(initial=None):
    """Native folder-selection dialog. Returns a path or None; raises if no dialog is available."""
    initial = initial or os.path.expanduser("~")
    try:
        import webview
        if getattr(webview, "windows", None):
            w = webview.windows[0]
            ft = getattr(getattr(webview, "FileDialog", None), "FOLDER", None)
            if ft is None:
                ft = getattr(webview, "FOLDER_DIALOG")
            r = w.create_file_dialog(ft, directory=initial)
            if isinstance(r, (list, tuple)):
                return r[0] if r else None
            return r or None
    except ImportError:
        pass
    except Exception:
        pass
    if sys.platform == "darwin":
        r = subprocess.run(["osascript", "-e", 'POSIX path of (choose folder with prompt '
                            '"Choose where CLAVAGEN saves results")'], capture_output=True, text=True)
        return r.stdout.strip() or None
    if sys.platform.startswith("win"):
        ps = ("Add-Type -AssemblyName System.Windows.Forms;"
              "$o=New-Object System.Windows.Forms.Form -Property @{TopMost=$true};"
              "$f=New-Object System.Windows.Forms.FolderBrowserDialog;"
              "$f.Description='Choose where CLAVAGEN saves results';"
              f"$f.SelectedPath='{initial}';"
              "if($f.ShowDialog($o) -eq 'OK'){ $f.SelectedPath }")
        r = subprocess.run(["powershell", "-STA", "-NoProfile", "-Command", ps], capture_output=True, text=True,
                           creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        return r.stdout.strip() or None
    for cmd in (["zenity", "--file-selection", "--directory", f"--filename={initial}/"],
                ["kdialog", "--getexistingdirectory", initial]):
        if shutil.which(cmd[0]):
            r = subprocess.run(cmd, capture_output=True, text=True)
            return r.stdout.strip() or None
    try:
        import tkinter
        from tkinter import filedialog
        root = tkinter.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        p = filedialog.askdirectory(initialdir=initial, title="Choose where CLAVAGEN saves results")
        root.destroy()
        return p or None
    except Exception:
        raise ValueError("no folder dialog available on this system")


def preview_annotation(data):
    from .online import auto_annotation, detect_gene_and_transcript
    from .seqio import SeqRecord
    ref = (data.get("reference") or "").strip()
    header = ref.splitlines()[0] if ref else ""
    ann = {}
    if (data.get("annotation") or "").strip():
        ann = json.loads(data["annotation"])
    rec = SeqRecord(id=header.lstrip(">").split(" ")[0] if header else "", seq="",
                    description=header.lstrip(">"))
    ctx = {"gene": (data.get("gene") or "").strip().upper() or ann.get("gene"), "transcript": ann.get("transcript")}
    m = re.fullmatch(r"(N[MR]_\d+\.\d+)", ref)
    if m:
        ctx["transcript"] = m.group(1)
    gene, tx = detect_gene_and_transcript(rec, ctx, resolve=True)
    if not gene and tx:
        from .online import fetch_refseq
        gene = fetch_refseq(tx).get("gene")
    if not gene:
        raise ValueError("Could not identify the gene: type it in 'Gene symbol', or put an NM_ accession / '(GENE)' "
                         "in the reference header.")
    a = auto_annotation(gene, tx)
    on = a.get("_online", {})
    kv = a.get("known_variants", [])
    return {"ok": True, "gene": gene, "transcript": a.get("transcript"), "inheritance": a.get("inheritance"),
            "lof": {True: "yes", False: "no"}.get(a.get("lof_mechanism"), "unknown"),
            "diseases": [f"{d['name']}" + (f" (MIM {d['mim']})" if d.get("mim") else "") for d in a.get("diseases", [])],
            "sources": on.get("sources", []) + [f"Derived: {x}" for x in a.get("_derived", [])],
            "errors": on.get("errors", []), "n_known": len(kv), "known": kv[:25]}


CLASS_ORDER = ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]


def _safe(s):
    return re.sub(r"[^\w.\-]+", "_", s)


def _multi_example():
    """Two genes in one run: HBB family (recessive) + SOD1 family (p.Asp91Ala, recessive form)."""
    from .examples import HBB_CDS, SOD1_CDS
    iu = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M"}

    def var(seq, pos, alt, het):
        s_ = list(seq)
        s_[pos - 1] = iu["".join(sorted(s_[pos - 1] + alt))] if het else alt
        return "\n".join("".join(s_)[i:i + 60] for i in range(0, len(s_), 60))
    wrap = lambda x: "\n".join(x[i:i + 60] for i in range(0, len(x), 60))  # noqa: E731
    refs = parse_reference_text(f">NM_000518.5 (HBB) hemoglobin subunit beta, coding sequence\n{wrap(HBB_CDS)}\n"
                                f">NM_000454.5 (SOD1) superoxide dismutase 1, coding sequence\n{wrap(SOD1_CDS)}\n")
    for r in refs:
        if r["gene"] == "SOD1":
            r["inheritance"] = "AD/AR"
    P = []

    def add(pid, fam, rel, aff, sex, age, ind, ph, seq, ref):
        P.append({"id": pid, "family": fam, "relationship": rel, "affected": aff, "sex": sex, "age": age,
                  "indication": ind, "phenotypes": ph, "evidence": "", "query": f">{pid}\n{seq}", "reference": ref,
                  "parentage_confirmed": False, "notes": "Demonstration data"})
    H = ("FAM-HBB-01", HBB_CDS, 20, "T", "NM_000518.5")
    for pid, rel, aff, sex, age, g in (("HBB_Father", "father", "unaffected", "male", "41", "het"),
                                       ("HBB_Mother", "mother", "unaffected", "female", "38", "het"),
                                       ("HBB_Child1", "proband", "affected", "female", "12", "hom"),
                                       ("HBB_Child2", "sibling", "affected", "male", "9", "hom"),
                                       ("HBB_Child3", "sibling", "unaffected", "female", "6", "het")):
        add(pid, H[0], rel, aff, sex, age, "Sickle cell family study", "haemolytic anaemia, vaso-occlusive pain crises"
            if aff == "affected" else "", var(H[1], H[2], H[3], g == "het"), "")
    S = ("FAM-SOD1-01", SOD1_CDS, 272, "C")
    for pid, rel, aff, sex, age, g in (("SOD1_Father", "father", "unaffected", "male", "68", "het"),
                                       ("SOD1_Mother", "mother", "unaffected", "female", "66", "het"),
                                       ("SOD1_Child1", "proband", "affected", "male", "42", "hom"),
                                       ("SOD1_Child2", "sibling", "unaffected", "female", "45", "het"),
                                       ("SOD1_Child3", "sibling", "unaffected", "male", "39", "het"),
                                       ("SOD1_Child4", "sibling", "unaffected", "female", "36", "wt"),
                                       ("SOD1_Child5", "sibling", "unaffected", "male", "33", "wt")):
        seq = wrap(S[1]) if g == "wt" else var(S[1], S[2], S[3], g == "het")
        add(pid, S[0], rel, aff, sex, age, "ALS family study", "progressive muscle weakness, fasciculations"
            if aff == "affected" else "", seq, "")
    return {"references": refs, "annotation": "", "known": "", "strict": "0", "patients": P}


def _family_example():
    """Trio + sibling (compound heterozygous + de novo) and a carrier couple, HBB."""
    base = _example_payload()
    from .examples import HBB_CDS
    iu = {"AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M"}

    def het(*ch):
        s = list(HBB_CDS)
        for pos, alt in ch:
            s[pos - 1] = iu["".join(sorted(s[pos - 1] + alt))]
        return "\n".join("".join(s)[i:i + 60] for i in range(0, len(s), 60))
    P = []

    def add(pid, fam, rel, aff, sex, age, ind, ph, seq, pc=False):
        P.append({"id": pid, "family": fam, "relationship": rel, "affected": aff, "sex": sex, "age": age,
                  "indication": ind, "phenotypes": ph, "evidence": "", "query": f">{pid}\n{seq}",
                  "parentage_confirmed": pc, "notes": "Demonstration data"})
    add("Child", "FAM1", "proband", "affected", "female", "4", "Severe anaemia, painful crises",
        "haemolytic anaemia, vaso-occlusive pain crises", het((20, "T"), (118, "T"), (250, "A")), True)
    add("Mother", "FAM1", "mother", "unaffected", "female", "33", "Parental testing", "", het((20, "T")))
    add("Father", "FAM1", "father", "unaffected", "male", "35", "Parental testing", "", het((118, "T")))
    add("Brother", "FAM1", "sibling", "unaffected", "male", "8", "Sibling testing", "", het((20, "T")))
    add("Woman_A", "COUPLE2", "proband", "unaffected", "female", "29", "Preconception carrier screening", "",
        het((20, "T")))
    add("Man_B", "COUPLE2", "partner", "unaffected", "male", "31", "Preconception carrier screening", "",
        het((19, "A")))
    base["patients"] = P
    return base


def _example_payload():
    import tempfile
    from .examples import write_examples
    from .seqio import parse_fasta_text
    d = tempfile.mkdtemp(prefix="acmg_example_")
    p = write_examples(d)
    rd = lambda k: open(p[k]).read()  # noqa: E731
    recs = parse_fasta_text(rd("dna_queries")) + parse_fasta_text(rd("protein_query"))
    info = {
        "Patient1_HBB": ("female", "24", "Anaemia work-up; family history of sickle cell disease",
                         "microcytic anaemia, vaso-occlusive pain crises", "c.20A>T=PS3,PS4"),
        "Patient2_HBB": ("male", "3", "Transfusion-dependent anaemia", "microcytic anaemia, hepatosplenomegaly", ""),
        "Patient3_HBB_amplicon": ("female", "31", "Antenatal carrier screening (partial amplicon)",
                                  "microcytic hypochromic anaemia", ""),
        "Patient4_HBB": ("male", "45", "Incidental finding follow-up", "", ""),
        "Patient5_HBB_protein": ("unknown", "", "Haemoglobin variant on mass spectrometry", "haemolytic anaemia", ""),
    }
    pts = []
    for rid, desc, seq in recs:
        sx, age, ind, ph, ev = info.get(rid, ("", "", "", "", ""))
        wrapped = "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
        pts.append({"id": rid.replace("_HBB", "").replace("_protein", "_prot").replace("_amplicon", "_amp"),
                    "sex": sx, "age": age, "indication": ind, "phenotypes": ph, "evidence": ev,
                    "query": f">{rid} {desc}\n{wrapped}", "notes": "Demonstration data"})
    return {"reference": rd("reference_dna"), "annotation": rd("annotation"), "gene": "", "cds": "", "strict": "0",
            "known": "", "patients": pts}


def _run_summary(folder, name):
    """Summary of a saved run for the Results / History views (multi-patient aware)."""
    from .cohort import patient_summary
    base = f"/runs/{name}/"
    rs = os.path.join(folder, "run_summary.json")
    patients = []
    if os.path.isfile(rs):
        with open(rs) as fh:
            meta = json.load(fh)
        date = meta["date"].replace("T", " ")
        for e in meta["patients"]:
            pbase = base + e["folder"] + "/"
            if not e.get("ok"):
                patients.append({"id": e["id"], "ok": False, "error": e.get("error", ""), "base": pbase,
                                 "samples": [], "n_variants": 0, "highest": None, "patient": {}, "phenotypes": []})
                continue
            ps = patient_summary(folder, e["folder"])
            for smp in ps["samples"]:
                smp["map"], smp["aln"] = pbase + smp["map"], pbase + smp["aln"]
            ps.update(id=e["id"], ok=True, base=pbase, report=pbase + "report.html")
            patients.append(ps)
        cohort = base + "cohort_summary.html"
        families = meta.get("families", [])
        for f in families:
            f["pedigree_url"] = base + f.get("pedigree", "")
    else:  # single-patient run from V3/V4
        ps = patient_summary(os.path.dirname(folder), os.path.basename(folder))
        for smp in ps["samples"]:
            smp["map"], smp["aln"] = base + smp["map"], base + smp["aln"]
        ps.update(id=(ps["patient"] or {}).get("id") or "Analysis", ok=True, base=base, report=base + "report.html")
        patients = [ps]
        date = ps["date"].replace("T", " ")
        cohort = None
        families = []
    allc = [p["highest"] for p in patients if p.get("highest")]
    highest = next((c for c in CLASS_ORDER if c in allc), None)
    ids = [p["id"] for p in patients]
    title = (ids[0] if len(ids) == 1 else f"{len(ids)} patients") + f" - {date}"
    return {"ok": True, "name": name, "base": base, "folder": os.path.abspath(folder), "date": date,
            "patients": patients, "cohort": cohort, "n_patients": len(patients), "families": families,
            "n_variants": sum(p["n_variants"] for p in patients), "highest": highest, "title": title,
            "patient_ids": ", ".join(ids)}


def parse_reference_text(text, label="pasted"):
    """FASTA/raw text or RefSeq accession(s) -> reference rows for the window."""
    from .online import detect_gene_and_transcript
    from .seqio import load_sequences
    text = (text or "").strip()
    if not text:
        return []
    lines = [x.strip() for x in text.splitlines() if x.strip()]
    if all(re.fullmatch(r"N[MR]_\d+(\.\d+)?", x) for x in lines):
        return [{"id": x, "accession": x, "description": "RefSeq accession (downloaded at run time)", "type": "dna",
                 "length": None, "gene": "", "transcript": x, "inheritance": "", "cds": "", "source": "accession"}
                for x in lines]
    import tempfile
    fd, tmp = tempfile.mkstemp(suffix=".fa")
    with os.fdopen(fd, "w") as fh:
        fh.write(text if text.startswith(">") else f">{label}\n{text}")
    rows = []
    for r in load_sequences([tmp], "ref"):
        gene, tx = detect_gene_and_transcript(r, {}, resolve=False)
        rows.append({"id": r.id, "description": r.description, "type": r.mol_type, "length": len(r.seq),
                     "gene": gene or "", "transcript": tx or "", "inheritance": "", "cds": "",
                     "fasta": f">{r.id} {r.description}".rstrip() + "\n" + r.seq, "source": "sequence"})
    os.remove(tmp)
    return rows


def _references_from_payload(data):
    """Reference list + per-reference annotation from the window payload (V9), or legacy text."""
    rows = data.get("references")
    if not rows:
        txt = (data.get("reference") or "").strip()
        if not txt:
            raise ValueError("Please add at least one reference (section 1).")
        rows = parse_reference_text(txt)
    ids = [r["id"] for r in rows]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        raise ValueError(f"Two references have the same ID ({', '.join(sorted(dup))}) - remove one or rename it.")
    refs, per = [], {}
    for r in rows:
        refs.append({"accession": r["accession"]} if r.get("accession") else {"fasta": r["fasta"], "id": r["id"]})
        st = {}
        if (r.get("gene") or "").strip():
            st["gene"] = r["gene"].strip().upper()
        if (r.get("transcript") or "").strip():
            st["transcript"] = r["transcript"].strip()
        if (r.get("inheritance") or "").strip():
            st["inheritance"] = r["inheritance"].strip()
        if (r.get("cds") or "").strip():
            a_, b_ = re.split(r"\s*(?:-|\.\.)\s*", r["cds"].strip())
            st["cds"] = [int(a_), int(b_)]
        if st:
            per[r["id"]] = st
    return refs, per


def run_payload(data, out_base):
    from .cohort import run_cohort
    from .online import parse_known_variants_text
    ann = {}
    if (data.get("annotation") or "").strip():
        try:
            ann = json.loads(data["annotation"])
        except json.JSONDecodeError as e:
            raise ValueError(f"Annotation is not valid JSON: {e}")
    refs, per = _references_from_payload(data)
    if per:
        ann = dict(ann)
        merged = dict(ann.get("references") or {})
        for rid, st in per.items():
            merged[rid] = dict(merged.get(rid) or {}, **st)   # settings in the reference list win
        ann["references"] = merged
    if (data.get("rules") or "").strip():
        import tempfile
        from .cli import _apply_rules
        try:
            json.loads(data["rules"])
        except json.JSONDecodeError as e:
            raise ValueError(f"Gene rule profile is not valid JSON: {e}")
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            fh.write(data["rules"])
        ann = _apply_rules(ann, fh.name)
        os.unlink(fh.name)
    known = parse_known_variants_text(data.get("known", ""), default_source="Manual entry")
    summ = run_cohort(data.get("patients") or [], refs, out_base, annotation=ann or None,
                      gene=(data.get("gene") or "").strip() or None, cds=(data.get("cds") or "").strip() or None,
                      strict=str(data.get("strict")) == "1", online=str(data.get("online", "1")) == "1",
                      known=known, use_clinvar=str(data.get("use_clinvar", "1")) == "1",
                      clinvar_min_stars=int(data.get("clinvar_min_stars") or 0))
    return _run_summary(os.path.join(out_base, summ["name"]), summ["name"])


def _open_path(p):
    if sys.platform.startswith("win"):
        os.startfile(p)  # noqa
    elif sys.platform == "darwin":
        subprocess.Popen(["open", p])
    else:
        subprocess.Popen(["xdg-open", p])


def make_handler(state, base_url):
    def cur_root():
        return os.path.realpath(state["results_dir"])

    def run_dir(name):
        root = cur_root()
        d = os.path.realpath(os.path.join(root, name or ""))
        if not d.startswith(root + os.sep) or not os.path.isdir(d):
            raise ValueError("Unknown run")
        return d

    class H(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _send(self, code, body, ctype="text/html; charset=utf-8"):
            if isinstance(body, str):
                body = body.encode()
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def _json(self, obj):
            return self._send(200, json.dumps(obj, default=str), "application/json")

        def do_GET(self):
            u = urlparse(self.path)
            path = unquote(u.path)
            try:
                if path == "/favicon.ico":
                    return self._send(204, b"", "image/x-icon")
                if path in ("/", "/index.html"):
                    return self._send(200, PAGE.replace("__VERSION__", __version__))
                if path == "/example":
                    q = u.query or ""
                    d = _multi_example() if "multi=1" in q else (_family_example() if "family=1" in q else _example_payload())
                    if "references" not in d:
                        d["references"] = parse_reference_text(d.get("reference", ""))
                    return self._json(d)
                if path == "/settings":
                    return self._json({"results_dir": state["results_dir"]})
                if path == "/history":
                    root = cur_root()
                    os.makedirs(root, exist_ok=True)
                    runs = []
                    for name in sorted(os.listdir(root), reverse=True):
                        d = os.path.join(root, name)
                        if os.path.isfile(os.path.join(d, "run_summary.json")) or \
                                os.path.isfile(os.path.join(d, "clinical_summary.json")):
                            try:
                                runs.append(_run_summary(d, name))
                            except Exception:
                                pass
                    return self._json({"folder": root, "runs": runs})
                if path == "/files":
                    d = run_dir(parse_qs(u.query).get("run", [""])[0])
                    files = []
                    for dp, _, fs in os.walk(d):
                        for f in sorted(fs):
                            full = os.path.join(dp, f)
                            files.append({"name": os.path.relpath(full, d).replace(os.sep, "/"),
                                          "size": f"{os.path.getsize(full) / 1024:.1f} KB"})
                    return self._json({"folder": d, "files": files})
                if path.startswith("/runs/"):
                    root = cur_root()
                    full = os.path.realpath(os.path.join(root, path[len("/runs/"):]))
                    if not full.startswith(root + os.sep) or not os.path.isfile(full):
                        return self._send(404, "Not found", "text/plain")
                    ctype = mimetypes.guess_type(full)[0] or "text/plain"
                    if full.endswith((".md", ".tsv", ".vcf", ".txt", ".fa")):
                        ctype = "text/plain; charset=utf-8"
                    if full.endswith(".svg"):
                        ctype = "image/svg+xml"
                    with open(full, "rb") as fh:
                        return self._send(200, fh.read(), ctype)
            except ValueError as e:
                return self._send(404, str(e), "text/plain")
            return self._send(404, "Not found", "text/plain")

        def do_POST(self):
            path = urlparse(self.path).path
            try:
                n = int(self.headers.get("Content-Length", 0))
                data = json.loads(self.rfile.read(n) or b"{}")
                if path == "/run":
                    if (data.get("results_dir") or "").strip() and \
                            os.path.abspath(os.path.expanduser(data["results_dir"].strip())) != state["results_dir"]:
                        state["results_dir"] = set_results_dir(data["results_dir"])
                    os.makedirs(state["results_dir"], exist_ok=True)
                    return self._json(run_payload(data, state["results_dir"]))
                if path == "/settings":
                    state["results_dir"] = set_results_dir(data.get("results_dir"))
                    return self._json({"ok": True, "results_dir": state["results_dir"]})
                if path == "/choosefolder":
                    p = choose_folder(state["results_dir"])
                    if p:
                        state["results_dir"] = set_results_dir(p)
                    return self._json({"ok": True, "path": state["results_dir"] if p else None})
                if path == "/checkconn":
                    from .online import check_connection
                    return self._json({"ok": True, "results": check_connection()})
                if path == "/annotate":
                    return self._json(preview_annotation(data))
                if path == "/parsereferences":
                    return self._json({"ok": True, "references": parse_reference_text(data.get("text", ""),
                                                                                       data.get("name") or "pasted")})
                if path == "/library":
                    lib = load_library()
                    act = data.get("action")
                    if act == "add":
                        for r in data.get("references", []):
                            lib = [x for x in lib if x["id"] != r["id"]] + [r]
                        save_library(lib)
                    elif act == "remove":
                        lib = [x for x in lib if x["id"] != data.get("id")]
                        save_library(lib)
                    return self._json({"ok": True, "library": lib})
                if path == "/parsepatients":
                    from .cohort import read_patient_table
                    pts = read_patient_table(data.get("text", ""))
                    for p in pts:
                        p["phenotypes"] = ", ".join(p["phenotypes"])
                        p["evidence"] = "; ".join(p["evidence"])
                    return self._json({"ok": True, "patients": pts})
                if path == "/savepng":
                    d = run_dir(data.get("run"))
                    target = os.path.realpath(os.path.join(d, data.get("file", "")))
                    if not target.endswith(".png") or not target.startswith(d + os.sep):
                        raise ValueError("bad file name")
                    raw = base64.b64decode(data["data"].split(",", 1)[1])
                    with open(target, "wb") as fh:
                        fh.write(raw)
                    return self._json({"ok": True})
                if path == "/open":
                    d = run_dir(data.get("run"))
                    if data.get("what") == "folder":
                        _open_path(d)
                    elif data.get("what") == "file":
                        full = os.path.realpath(os.path.join(d, data.get("name", "")))
                        if not full.startswith(d + os.sep) or not os.path.isfile(full):
                            raise ValueError("Unknown file")
                        _open_path(full)
                    else:
                        relp = (data.get("path") or "report.html").lstrip("/")
                        full = os.path.realpath(os.path.join(d, relp))
                        if not full.startswith(d + os.sep) or not os.path.isfile(full):
                            raise ValueError("report not found")
                        webbrowser.open(f"{base_url}runs/{os.path.basename(d)}/{relp}")
                    return self._json({"ok": True})
            except Exception as e:  # every error is reported back to the window
                return self._json({"ok": False, "error": str(e)})
            return self._send(404, "Not found", "text/plain")
    return H


def _free_port(preferred=8765):
    for p in [preferred] + list(range(preferred + 1, preferred + 50)):
        with socket.socket() as s:
            try:
                s.bind(("127.0.0.1", p))
                return p
            except OSError:
                continue
    return 0


# ------------------------------------------------------------------ windows
def _find_app_browser():
    """Chromium-family browser that supports --app (chrome-less application window)."""
    cands = []
    if sys.platform.startswith("win"):
        for env in ("PROGRAMFILES(X86)", "PROGRAMFILES", "LOCALAPPDATA"):
            b = os.environ.get(env)
            if b:
                cands += [os.path.join(b, "Microsoft", "Edge", "Application", "msedge.exe"),
                          os.path.join(b, "Google", "Chrome", "Application", "chrome.exe"),
                          os.path.join(b, "BraveSoftware", "Brave-Browser", "Application", "brave.exe")]
    elif sys.platform == "darwin":
        for app, exe in (("Google Chrome", "Google Chrome"), ("Microsoft Edge", "Microsoft Edge"),
                         ("Brave Browser", "Brave Browser"), ("Chromium", "Chromium")):
            for base in ("/Applications", os.path.expanduser("~/Applications")):
                cands.append(f"{base}/{app}.app/Contents/MacOS/{exe}")
    else:
        for n in ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser", "microsoft-edge", "brave-browser"):
            w = shutil.which(n)
            if w:
                cands.append(w)
    return next((c for c in cands if os.path.isfile(c)), None)


def _native_window(url):
    import webview  # pywebview
    webview.create_window(f"{APP_TITLE} {__version__}", url, width=1320, height=900, min_size=(960, 640))
    webview.start()


def _app_mode_window(url, exe):
    profile = app_path("window_profile")
    os.makedirs(profile, exist_ok=True)
    proc = subprocess.Popen([exe, f"--app={url}", f"--user-data-dir={profile}", "--no-first-run",
                             "--no-default-browser-check", "--window-size=1320,900"])
    proc.wait()


def main(out_base=None, port=None, open_browser=True, window="auto"):
    """out_base: explicit results folder (else the last folder chosen in the window, else
    ~/Documents/CLAVAGEN Results)."""
    out_base = os.path.abspath(out_base) if out_base else (load_settings().get("results_dir") or default_results_dir())
    os.makedirs(out_base, exist_ok=True)
    state = {"results_dir": out_base}
    port = port or _free_port()
    url_holder = {}
    srv = ThreadingHTTPServer(("127.0.0.1", port), None)
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    url_holder["url"] = url
    srv.RequestHandlerClass = make_handler(state, url)
    print(f"{APP_TITLE} {__version__} - results are saved in: {out_base}")
    if not open_browser:
        print(f"Serving at {url} (Ctrl+C to stop)")
        try:
            srv.serve_forever()
        except KeyboardInterrupt:
            pass
        return
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    try:
        if window in ("auto", "native"):
            try:
                _native_window(url)
                return
            except Exception as e:
                if window == "native":
                    raise
                print(f"(native window unavailable: {e.__class__.__name__}: {e})")
        exe = _find_app_browser() if window in ("auto", "app") else None
        if exe:
            print(f"Opening application window ({os.path.basename(exe)} app mode). Close the window to quit.")
            _app_mode_window(url, exe)
            return
        print(f"Opening {url} in the default browser. Press Ctrl+C here to quit.")
        webbrowser.open(url)
        while t.is_alive():
            t.join(1)
    except KeyboardInterrupt:
        pass
    finally:
        srv.shutdown()
        srv.server_close()


if __name__ == "__main__":
    main()
