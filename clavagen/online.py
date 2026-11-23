"""Automatic online annotation (standard library only).

Sources (all public, no API key):
* NCBI ClinVar (E-utilities)   - known pathogenic / likely pathogenic variants, MANE/RefSeq transcript
* NCBI Nucleotide (E-utilities) - RefSeq transcript sequence, CDS and exon structure
* UniProt REST                 - protein domains, function, disease associations (MIM)
* gnomAD GraphQL               - gene constraint (LOEUF, missense Z)
* Ensembl VEP REST             - per-variant gnomAD allele frequencies, REVEL, SpliceAI, AlphaMissense

Every response is cached in ~/.clavagen/cache (default 30 days) so analyses are
reproducible and work offline after the first lookup. All failures are non-fatal:
the run continues with whatever annotation is available, and the reasons are reported.
Data sent: gene symbols, transcript accessions and variant HGVS strings only - never
patient identifiers or full patient sequences.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter

from . import __version__
from .genetic_code import AA1

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
GNOMAD = "https://gnomad.broadinstitute.org/api"
VEP = "https://rest.ensembl.org/vep/human/hgvs"
from .paths import app_path  # noqa: E402
CACHE_DIR = app_path("cache")
CACHE_DAYS = 30
UA = f"CLAVAGEN/{__version__} (clinical decision support; python-urllib)"

_last_ncbi = [0.0]
_dead_hosts = {}  # host -> reason; hosts that could not be reached in this session
_SSL = [None]


def ssl_context():
    """TLS context that also works with python.org macOS builds (which ship without CA
    certificates): certifi if installed, else the system bundle."""
    if _SSL[0] is not None:
        return _SSL[0]
    import ssl
    import sys
    ctx = None
    try:
        import certifi
        ctx = ssl.create_default_context(cafile=certifi.where())
    except Exception:
        ctx = ssl.create_default_context()
        for bundle in ("/etc/ssl/cert.pem", "/usr/local/etc/openssl/cert.pem", "/opt/homebrew/etc/openssl@3/cert.pem",
                       "/etc/ssl/certs/ca-certificates.crt", "/etc/pki/tls/certs/ca-bundle.crt"):
            if os.path.isfile(bundle):
                try:
                    ctx.load_verify_locations(bundle)
                except Exception:
                    pass
        if sys.platform.startswith("win"):
            try:
                ctx.load_default_certs()
            except Exception:
                pass
    _SSL[0] = ctx
    return ctx


def explain_error(e):
    """Turn a network exception into an actionable message for clinicians."""
    import socket
    import ssl
    txt = str(e)
    reason = getattr(e, "reason", None)
    if isinstance(e, ssl.SSLError) or isinstance(reason, ssl.SSLError) or "CERTIFICATE_VERIFY_FAILED" in txt:
        return ("secure connection failed (SSL certificate check). On macOS open Applications > Python 3.x > "
                "'Install Certificates.command', or start CLAVAGEN with its start file (it installs certificates).")
    if isinstance(e, urllib.error.HTTPError):
        return f"server answered HTTP {e.code} ({e.reason})"
    if isinstance(reason, socket.gaierror) or "Name or service not known" in txt or "nodename nor servname" in txt:
        return "no internet connection (server name could not be resolved)"
    if isinstance(e, (socket.timeout, TimeoutError)) or "timed out" in txt:
        return "connection timed out (slow network, firewall or proxy)"
    if "403" in txt or "Forbidden" in txt or "Tunnel connection failed" in txt:
        return "blocked by a firewall/proxy (institutional network?)"
    return txt


class OnlineError(Exception):
    pass


# ------------------------------------------------------------------ http
def _cache_path(key):
    return os.path.join(CACHE_DIR, hashlib.sha256(key.encode()).hexdigest()[:32] + ".json")


def http(url, data=None, headers=None, timeout=30, use_cache=True, kind="json"):
    """GET/POST with on-disk cache. Returns parsed JSON (kind='json') or text."""
    key = url + "|" + (json.dumps(data, sort_keys=True) if data is not None else "")
    cp = _cache_path(key)
    if use_cache and os.path.isfile(cp) and (time.time() - os.path.getmtime(cp)) < CACHE_DAYS * 86400:
        with open(cp) as fh:
            return json.load(fh)["body"]
    host = urllib.parse.urlparse(url).netloc
    if host in _dead_hosts:
        raise OnlineError(f"{host}: {_dead_hosts[host]}")
    if url.startswith(EUTILS):  # NCBI asks for <= 3 requests/second without an API key
        wait = 0.35 - (time.time() - _last_ncbi[0])
        if wait > 0:
            time.sleep(wait)
        _last_ncbi[0] = time.time()
    hdr = {"User-Agent": UA, "Accept": "application/json" if kind == "json" else "*/*"}
    hdr.update(headers or {})
    body = None
    if data is not None:
        body = json.dumps(data).encode()
        hdr["Content-Type"] = "application/json"
    last = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=body, headers=hdr)
            with urllib.request.urlopen(req, timeout=timeout, context=ssl_context()) as r:
                raw = r.read().decode("utf-8", errors="replace")
            out = json.loads(raw) if kind == "json" else raw
            if use_cache:
                os.makedirs(CACHE_DIR, exist_ok=True)
                with open(cp, "w") as fh:
                    json.dump({"url": url, "retrieved": _dt.datetime.now().isoformat(timespec="seconds"),
                               "body": out}, fh)
            return out
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(1.5 * (attempt + 1))
                continue
            break
        except (urllib.error.URLError, OSError) as e:  # network down, DNS, SSL, proxy refusal, timeout
            last = e
            msg = explain_error(e)
            if "SSL" in msg or "internet" in msg or attempt == 1:
                _dead_hosts[host] = msg
                break
            time.sleep(1.0)
        except Exception as e:  # bad JSON etc.
            last = e
            break
    raise OnlineError(f"{host}: {explain_error(last) if last else 'unknown error'}")


def _q(params):
    return urllib.parse.urlencode(params)


# ------------------------------------------------------------------ ClinVar
_TITLE = re.compile(r"^(?P<tx>N[MR]_\d+(?:\.\d+)?)\((?P<gene>[A-Za-z0-9\-]+)\):(?P<c>c\.[^\s]+)(?:\s+\((?P<p>p\.[^)]+)\))?")


REVIEW_STARS = [("practice guideline", 4), ("expert panel", 3), ("multiple submitters", 2),
                ("conflicting", 1), ("single submitter", 1), ("criteria provided", 1),
                ("no assertion", 0), ("no classification", 0)]


def review_stars(status):
    st = (status or "").lower()
    for key, n in REVIEW_STARS:
        if key in st:
            return n
    return 0


def _classification(doc):
    for k in ("germline_classification", "clinical_significance"):
        v = doc.get(k)
        if isinstance(v, dict) and v.get("description"):
            return v.get("description"), v.get("review_status", "")
    return "", ""


def _traits(doc):
    names = []
    for k in ("germline_classification", "clinical_significance"):
        v = doc.get(k)
        if isinstance(v, dict):
            for t in v.get("trait_set") or []:
                n = t.get("trait_name") if isinstance(t, dict) else None
                if n and n not in names:
                    names.append(n)
    for t in doc.get("trait_set") or []:
        n = t.get("trait_name") if isinstance(t, dict) else None
        if n and n not in names:
            names.append(n)
    return [n for n in names if n.lower() not in ("not provided", "not specified")]


def parse_clinvar_summary(result, gene):
    """Convert an esummary 'result' block to known_variants entries (P/LP only)."""
    out = []
    for uid in result.get("uids", []):
        doc = result.get(uid) or {}
        title = doc.get("title", "")
        m = _TITLE.match(title)
        if not m or m.group("gene").upper() != gene.upper():
            continue
        cls, rev = _classification(doc)
        cl = cls.lower()
        if not (cl.startswith("pathogenic") or cl.startswith("likely pathogenic")):
            continue
        if "conflicting" in cl:
            continue
        p = m.group("p") or doc.get("protein_change") or ""
        if p and not p.startswith("p."):
            p = _one_letter_to_hgvs(p)
        out.append({"hgvs_c": m.group("c"), "hgvs_p": p, "transcript": m.group("tx"),
                    "classification": cls, "review_status": rev, "stars": review_stars(rev),
                    "conditions": _traits(doc),
                    "source": f"ClinVar VCV{int(uid):09d}" if str(uid).isdigit() else "ClinVar"})
    return out


def _one_letter_to_hgvs(pc):
    # ClinVar 'protein_change' like "E7V" (may be comma-separated list)
    pc = pc.split(",")[0].strip()
    m = re.match(r"^([A-Z])(\d+)([A-Z*])$", pc)
    if not m:
        return ""
    from .genetic_code import AA3
    return f"p.{AA3.get(m.group(1))}{m.group(2)}{AA3.get(m.group(3), 'Ter')}"


def fetch_clinvar(gene, max_records=3000, log=None):
    log = log or (lambda *_: None)
    term = f'{gene}[gene] AND ("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])'
    r = http(EUTILS + "esearch.fcgi?" + _q({"db": "clinvar", "term": term, "retmax": max_records,
                                             "retmode": "json"}))
    ids = r.get("esearchresult", {}).get("idlist", [])
    if not ids:  # fall back to all records for the gene, filtered locally
        r = http(EUTILS + "esearch.fcgi?" + _q({"db": "clinvar", "term": f"{gene}[gene]", "retmax": max_records,
                                                 "retmode": "json"}))
        ids = r.get("esearchresult", {}).get("idlist", [])
    log(f"ClinVar: {len(ids)} records for {gene}")
    out = []
    for i in range(0, len(ids), 200):
        s = http(EUTILS + "esummary.fcgi?" + _q({"db": "clinvar", "id": ",".join(ids[i:i + 200]),
                                                  "retmode": "json"}))
        out += parse_clinvar_summary(s.get("result", {}), gene)
    # de-duplicate on c.
    seen, uniq = set(), []
    for kv in out:
        if kv["hgvs_c"] not in seen:
            seen.add(kv["hgvs_c"])
            uniq.append(kv)
    return uniq


# ------------------------------------------------------------------ RefSeq
def parse_genbank(text):
    """Minimal GenBank parser: sequence, gene, CDS span, protein_id, exon spans."""
    acc = re.search(r"^VERSION\s+(\S+)", text, re.M)
    seq_part = text.split("\nORIGIN", 1)[1] if "\nORIGIN" in text else ""
    seq = re.sub(r"[^acgtnACGTN]", "", seq_part.split("//")[0]).upper()
    feats = text.split("\nFEATURES", 1)[1].split("\nORIGIN", 1)[0] if "\nFEATURES" in text else ""
    cds = re.search(r"^\s{5}CDS\s+(?:join\()?<?(\d+)\.\.>?(\d+)", feats, re.M)
    exons = [(int(a), int(b)) for a, b in re.findall(r"^\s{5}exon\s+<?(\d+)\.\.>?(\d+)", feats, re.M)]
    gene = re.search(r'/gene="([^"]+)"', feats)
    prot = re.search(r'/protein_id="([^"]+)"', feats)
    out = {"transcript": acc.group(1) if acc else None, "seq": seq, "gene": gene.group(1) if gene else None,
           "protein_id": prot.group(1) if prot else None}
    if cds:
        cs, ce = int(cds.group(1)), int(cds.group(2))
        out["cds"] = [cs, ce]
        # exon junctions in c. coordinates (last base of each exon except the last)
        out["exon_boundaries_c"] = [e - cs + 1 for _, e in exons[:-1] if cs <= e < ce]
        out["n_exons"] = len(exons)
    return out


def fetch_protein_record(accession):
    """NP_ protein record -> {gene, transcript (coded_by NM_), protein_id, seq}."""
    txt = http(EUTILS + "efetch.fcgi?" + _q({"db": "protein", "id": accession, "rettype": "gp",
                                              "retmode": "text"}), kind="text")
    if "LOCUS" not in txt:
        raise OnlineError(f"RefSeq protein {accession} not found")
    gene = re.search(r'/gene="([^"]+)"', txt)
    coded = re.search(r'/coded_by="(?:complement\()?(N[MR]_\d+\.\d+)', txt)
    ver = re.search(r"^VERSION\s+(\S+)", txt, re.M)
    seq_part = txt.split("\nORIGIN", 1)[1].split("//")[0] if "\nORIGIN" in txt else ""
    return {"gene": gene.group(1) if gene else None, "transcript": coded.group(1) if coded else None,
            "protein_id": ver.group(1) if ver else accession,
            "seq": re.sub(r"[^A-Za-z]", "", seq_part).upper()}


def fetch_refseq(accession):
    txt = http(EUTILS + "efetch.fcgi?" + _q({"db": "nuccore", "id": accession, "rettype": "gb",
                                              "retmode": "text"}), kind="text")
    if "LOCUS" not in txt:
        raise OnlineError(f"RefSeq record {accession} not found")
    return parse_genbank(txt)


# ------------------------------------------------------------------ UniProt
def parse_uniprot(js):
    res = (js or {}).get("results") or []
    if not res:
        return {}
    e = res[0]
    out = {"uniprot": e.get("primaryAccession"), "domains": [], "diseases": [], "function": "",
           "critical_sites": []}
    for f in e.get("features", []):
        t = f.get("type", "")
        loc = f.get("location", {})
        s = (loc.get("start") or {}).get("value")
        en = (loc.get("end") or {}).get("value")
        if s is None or en is None:
            continue
        if t in ("Domain", "Zinc finger", "DNA binding", "Repeat", "Motif", "Region", "Transmembrane"):
            if t == "Region" and "disordered" in (f.get("description") or "").lower():
                continue
            out["domains"].append({"name": f.get("description") or t, "start": int(s), "end": int(en),
                                   "critical": False, "type": t, "source": "UniProt"})
        elif t in ("Active site", "Binding site", "Site", "Metal binding"):
            out["critical_sites"] += list(range(int(s), int(en) + 1))
    for c in e.get("comments", []):
        ct = c.get("commentType")
        if ct == "FUNCTION" and c.get("texts"):
            out["function"] = c["texts"][0].get("value", "")
        if ct == "DISEASE" and c.get("disease"):
            d = c["disease"]
            desc = d.get("description", "")
            low = desc.lower()
            modes = []
            for key, code in (("autosomal dominant", "AD"), ("autosomal recessive", "AR"),
                              ("x-linked dominant", "XLD"), ("x-linked recessive", "XLR")):
                if key in low:
                    modes.append(code)
            if not any(m.startswith("XL") for m in modes) and "x-linked" in low:
                modes.append("XL")
            inh = "/".join(modes)
            xref = d.get("diseaseCrossReference") or {}
            out["diseases"].append({"name": d.get("diseaseId", ""), "mim": xref.get("id", "") if xref.get("database") == "MIM" else "",
                                    "inheritance": inh, "summary": desc, "features": [], "source": "UniProt"})
    return out


def fetch_uniprot(gene):
    q = f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true"
    js = http(UNIPROT + "?" + _q({"query": q, "format": "json", "size": 1}))
    return parse_uniprot(js)


# ------------------------------------------------------------------ gnomAD
def fetch_gnomad_constraint(gene):
    q = ("query($g:String!){gene(gene_symbol:$g, reference_genome:GRCh38){symbol "
         "gnomad_constraint{oe_lof_upper mis_z}}}")
    js = http(GNOMAD, data={"query": q, "variables": {"g": gene}})
    g = ((js or {}).get("data") or {}).get("gene") or {}
    c = g.get("gnomad_constraint") or {}
    return {"loeuf": c.get("oe_lof_upper"), "mis_z": c.get("mis_z")}


# ------------------------------------------------------------------ VEP
def _walk(o):
    if isinstance(o, dict):
        for k, v in o.items():
            yield k, v
            yield from _walk(v)
    elif isinstance(o, list):
        for v in o:
            yield from _walk(v)


def parse_vep_record(rec, alt=None, transcript=None):
    """Extract gnomAD AF, ClinVar significance and predictor scores from one VEP result."""
    out = {}
    afs = {}
    for cv in rec.get("colocated_variants", []) or []:
        freqs = cv.get("frequencies") or {}
        alleles = [alt] if alt and alt in freqs else list(freqs)
        for al in alleles:
            for k, v in (freqs.get(al) or {}).items():
                if k in ("gnomade", "gnomadg", "gnomad") and isinstance(v, (int, float)):
                    afs[k] = max(afs.get(k, 0), float(v))
        if cv.get("clin_sig"):
            out["clinvar_vep"] = ",".join(cv["clin_sig"])
        if cv.get("id", "").startswith("rs"):
            out["rsid"] = cv["id"]
    if afs:
        out["af"] = max(afs.values())
        out["af_detail"] = afs
    tcs = rec.get("transcript_consequences", []) or []
    if transcript:
        base = transcript.split(".")[0]
        sel = [t for t in tcs if str(t.get("transcript_id", "")).split(".")[0] == base]
        tcs = sel or tcs
    revel = spl = am = None
    for k, v in _walk(tcs):
        lk = k.lower()
        if lk in ("revel", "revel_score") and isinstance(v, (int, float, str)):
            try:
                revel = max(revel or 0, float(v))
            except ValueError:
                pass
        elif lk == "am_pathogenicity":
            try:
                am = max(am or 0, float(v))
            except (TypeError, ValueError):
                pass
        elif lk.startswith("spliceai_pred_ds") or lk in ("ds_ag", "ds_al", "ds_dg", "ds_dl"):
            try:
                spl = max(spl or 0, float(v))
            except (TypeError, ValueError):
                pass
    if revel is not None:
        out["REVEL"] = revel
    if spl is not None:
        out["SpliceAI"] = spl
    if am is not None:
        out["AlphaMissense"] = am
    return out


def vep_lookup(hgvs_list, transcript=None):
    """hgvs_list: ['NM_000518.5:c.20A>T', ...] -> {hgvs: parsed dict}."""
    res = {}
    for i in range(0, len(hgvs_list), 50):
        chunk = hgvs_list[i:i + 50]
        js = http(VEP + "?" + _q({"refseq": 1, "REVEL": 1, "SpliceAI": 1, "AlphaMissense": 1,
                                  "content-type": "application/json"}),
                  data={"hgvs_notations": chunk}, timeout=60)
        for rec in js or []:
            inp = rec.get("input", "")
            alt = None
            m = re.search(r">([ACGT])$", inp)
            if m:
                alt = m.group(1)
            res[inp] = parse_vep_record(rec, alt, transcript)
    return res


# ------------------------------------------------------------------ build
def detect_gene_and_transcript(ref_record, ctx, resolve=True, log=None):
    """Gene symbol and transcript from: annotation/gene box, NM_ accession, '(GENE)', 'gene=GENE',
    'GN=GENE' (UniProt) in the header, or - online - an NP_ protein accession (via NCBI Protein)."""
    tx = ctx.get("transcript")
    gene = ctx.get("gene")
    text = f"{ref_record.id} {ref_record.description}"
    if not tx:
        m = re.search(r"(N[MR]_\d+\.\d+)", text)
        tx = m.group(1) if m else None
    npm = re.search(r"(NP_\d+\.\d+)", text)
    if npm and resolve and (not gene or not tx):
        try:
            pr = fetch_protein_record(npm.group(1))
            gene = gene or pr.get("gene")
            tx = tx or pr.get("transcript")
            if log:
                log(f"{npm.group(1)} -> gene {pr.get('gene')}, transcript {pr.get('transcript')} (NCBI Protein)")
        except Exception as e:
            if log:
                log(f"Could not resolve {npm.group(1)}: {e}")
    if not gene:
        m = (re.search(r"\(([A-Z0-9\-]{2,15})\)", text) or re.search(r"\bgene[=:]\s*([A-Za-z0-9\-]+)", text, re.I)
             or re.search(r"\bGN=([A-Za-z0-9\-]+)", text))
        if m:
            gene = m.group(1).upper()
        else:
            for tok in re.split(r"[\s_|,;]+", text):
                if re.fullmatch(r"[A-Z][A-Z0-9]{1,9}", tok) and not re.match(r"N[MRPG]$", tok):
                    from .knowledge import GENES
                    if tok in GENES:
                        gene = tok
                        break
    return gene, tx


def auto_annotation(gene, transcript=None, log=None, sources=("clinvar", "refseq", "uniprot", "gnomad")):
    """Build an annotation dict for a gene from public databases. Never raises."""
    log = log or (lambda *_: None)
    ann = {"gene": gene, "_online": {"retrieved": _dt.datetime.now().isoformat(timespec="seconds"),
                                     "sources": [], "errors": []}}
    meta = ann["_online"]
    known = []
    if "clinvar" in sources:
        try:
            known = fetch_clinvar(gene, log=log)
            ann["known_variants"] = known
            meta["sources"].append(f"ClinVar (E-utilities): {len(known)} P/LP variants")
            log(f"ClinVar: {len(known)} pathogenic/likely pathogenic variants")
        except Exception as e:
            meta["errors"].append(f"ClinVar: {e}")
    # transcript: given, else most common in ClinVar titles (usually MANE Select)
    if not transcript and known:
        transcript = Counter(k["transcript"] for k in known if k.get("transcript")).most_common(1)[0][0]
    if transcript:
        ann["transcript"] = transcript
        known_tx = [k for k in known if not k.get("transcript") or k["transcript"].split(".")[0] == transcript.split(".")[0]]
        if known:
            ann["known_variants"] = known_tx
    if "refseq" in sources and transcript:
        try:
            gb = fetch_refseq(transcript)
            ann["_refseq"] = gb
            ann["protein_id"] = gb.get("protein_id")
            if gb.get("exon_boundaries_c") is not None:
                ann["exon_boundaries_c"] = gb["exon_boundaries_c"]
            meta["sources"].append(f"NCBI RefSeq {transcript}: CDS {gb.get('cds')}, {gb.get('n_exons', '?')} exons")
        except Exception as e:
            meta["errors"].append(f"RefSeq: {e}")
    if "uniprot" in sources:
        try:
            up = fetch_uniprot(gene)
            if up:
                ann["domains"] = up["domains"]
                if up["diseases"]:
                    ann["diseases"] = up["diseases"]
                    inh = {m for d in up["diseases"] for m in d["inheritance"].split("/") if m}
                    if len(inh) == 1:
                        ann["inheritance"] = inh.pop()
                    elif inh:
                        ann["inheritance"] = "/".join(sorted(inh))
                ann["function"] = up.get("function", "")
                ann["_critical_sites"] = up["critical_sites"]
                meta["sources"].append(f"UniProt {up.get('uniprot')}: {len(up['domains'])} domains/regions, "
                                       f"{len(up['diseases'])} disease(s)")
        except Exception as e:
            meta["errors"].append(f"UniProt: {e}")
    if "gnomad" in sources:
        try:
            gc = fetch_gnomad_constraint(gene)
            ann["constraint"] = gc
            if gc.get("mis_z") is not None and gc["mis_z"] >= 3.09:
                ann["missense_constrained"] = True
            meta["sources"].append(f"gnomAD constraint: LOEUF {gc.get('loeuf')}, missense Z {gc.get('mis_z')}")
        except Exception as e:
            meta["errors"].append(f"gnomAD: {e}")
    derive_gene_evidence(ann)
    return ann


def derive_gene_evidence(ann):
    """Derive LOF mechanism and hotspots from the known P/LP variants and constraint."""
    known = ann.get("known_variants") or []
    trunc = [k for k in known if re.search(r"(Ter|\*|fs)", k.get("hgvs_p", "")) or
             re.search(r"c\.\d+[+-][12][ACGT]>", k.get("hgvs_c", ""))]
    loeuf = (ann.get("constraint") or {}).get("loeuf")
    notes = ann.setdefault("_derived", [])
    if "lof_mechanism" not in ann:
        if len(trunc) >= 3:
            ann["lof_mechanism"] = True
            ann["_lof_derived"] = "clinvar"     # remembered so leave-one-out can recount without the target
            notes.append(f"LOF mechanism inferred: {len(trunc)} P/LP truncating/canonical-splice variants in ClinVar.")
        elif loeuf is not None and loeuf < 0.35:
            ann["lof_mechanism"] = True
            notes.append(f"LOF mechanism inferred from gnomAD LOEUF {loeuf:.2f} (<0.35).")
    # hotspots: residues with >= 2 different P/LP missense changes
    pos = Counter()
    for k in known:
        m = re.match(r"^p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)?$", k.get("hgvs_p", ""))
        if m and m.group(3) not in ("Ter",) and AA1.get(m.group(1)):
            pos[int(m.group(2))] += 1
    hs = sorted(p for p, n in pos.items() if n >= 2)
    if hs:
        ann["_hotspots_derived"] = sorted(set(hs) - set(ann.get("hotspots", [])))
        ann["hotspots"] = sorted(set(ann.get("hotspots", [])) | set(hs))
        notes.append(f"{len(hs)} missense hotspot residue(s) (>=2 distinct P/LP missense in ClinVar).")
    return ann


def recalc_gene_evidence(ctx, fixed_hotspots=(), fixed_lof=None):
    """Leave-one-out (benchmarking): recount the ClinVar-derived gene-level evidence after the tested
    variant has been removed from ctx['known_variants'].

    Hotspots ('>=2 distinct P/LP missense at a residue') and the loss-of-function mechanism ('>=3 P/LP
    truncating variants') are derived from ClinVar before the exclusion, so the tested variant could
    count towards its own PM1/PVS1 evidence. Values from the built-in gene list or the user's annotation
    (``fixed_hotspots``, ``fixed_lof``) are kept as they are.
    """
    before_hs = set(ctx.get("hotspots") or [])
    before_lof = ctx.get("lof_mechanism")
    prior_derived = set(ctx.get("_hotspots_derived") or [])
    tmp = {"known_variants": ctx.get("known_variants") or [], "constraint": ctx.get("constraint") or {}}
    if fixed_lof is not None:
        tmp["lof_mechanism"] = fixed_lof
    derive_gene_evidence(tmp)
    keep = (before_hs - prior_derived) | set(fixed_hotspots or [])
    # removing a variant can only take hotspots away, never add ones that were not in effect before
    ctx["hotspots"] = sorted(keep | (set(tmp.get("_hotspots_derived") or []) & before_hs))
    ctx["_hotspots_derived"] = tmp.get("_hotspots_derived") or []
    if ctx.get("_lof_derived") == "clinvar" and fixed_lof is None:
        ctx["lof_mechanism"] = tmp.get("lof_mechanism")
        ctx["_lof_derived"] = tmp.get("_lof_derived")
    notes = []
    lost = sorted(before_hs - set(ctx["hotspots"]))
    if lost:
        notes.append(f"Leave-one-out: hotspot residue(s) {', '.join(map(str, lost))} no longer qualify without the "
                     "tested variant.")
    if bool(before_lof) and not ctx.get("lof_mechanism"):
        notes.append("Leave-one-out: the loss-of-function mechanism is no longer supported without the tested variant.")
    ctx.setdefault("_derived", [])
    ctx["_derived"] = list(ctx["_derived"]) + notes
    return notes


def merge_annotation(base, auto, manual_known=None):
    """auto < built-in knowledge < user annotation file; manual known variants appended."""
    out = dict(auto or {})
    for k, v in (base or {}).items():
        if k == "known_variants":
            continue
        if k == "references":
            out[k] = v
            continue
        out[k] = v
    kv = list((auto or {}).get("known_variants", [])) + list((base or {}).get("known_variants", []))
    kv += list(manual_known or [])
    seen, merged = set(), []
    for x in reversed(kv):  # later (manual/user) entries take precedence
        key = (x.get("hgvs_c") or "") + "|" + (x.get("hgvs_p") or "")
        if key not in seen:
            seen.add(key)
            merged.append(x)
    out["known_variants"] = list(reversed(merged))
    return out


def parse_known_variants_text(text, default_source="User list"):
    """Parse manual known-variant entries.

    One per line; separators: tab, comma, semicolon or '|'. Columns may be any of
    c. notation, p. notation, classification, source (order free; detected by pattern).
    Lines starting with '#' and header lines are ignored. Also accepts JSON lists.
    """
    text = (text or "").strip()
    if not text:
        return []
    if text.startswith("["):
        return [dict(x, source=x.get("source", default_source)) for x in json.loads(text)]
    out = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.lower().startswith(("hgvs", "variant", "c.\t")):
            continue
        parts = [p.strip() for p in re.split(r"[\t,;|]", line) if p.strip()]
        rec = {"hgvs_c": "", "hgvs_p": "", "classification": "Pathogenic", "source": default_source}
        rest = []
        for p in parts:
            pp = p.split(":", 1)[1] if re.match(r"^N[MP]_\d", p) and ":" in p else p
            if pp.startswith("c."):
                rec["hgvs_c"] = pp
            elif pp.startswith("p."):
                rec["hgvs_p"] = pp
            elif re.match(r"^(likely\s+)?pathogenic|^(likely\s+)?benign|^uncertain|^vus$|^p$|^lp$", pp, re.I):
                rec["classification"] = {"p": "Pathogenic", "lp": "Likely pathogenic", "vus": "Uncertain significance"}.get(
                    pp.lower(), pp[0].upper() + pp[1:])
            elif re.match(r"^[A-Z]\d+[A-Z*]$", pp):
                rec["hgvs_p"] = _one_letter_to_hgvs(pp)
            else:
                rest.append(p)
        if rest:
            rec["source"] = " ".join(rest)
        if rec["hgvs_c"] or rec["hgvs_p"]:
            out.append(rec)
    return out


def check_connection():
    """Quick reachability test of each service (no cache). Returns list of (service, ok, message)."""
    tests = [
        ("NCBI (ClinVar / RefSeq)", EUTILS + "einfo.fcgi?db=clinvar&retmode=json", None),
        ("UniProt", UNIPROT + "?query=gene_exact:SOD1+AND+organism_id:9606&format=json&size=1&fields=accession", None),
        ("gnomAD", GNOMAD, {"query": "{__typename}"}),
        ("Ensembl VEP", "https://rest.ensembl.org/info/ping?content-type=application/json", None),
    ]
    out = []
    for name, url, data in tests:
        host = urllib.parse.urlparse(url).netloc
        _dead_hosts.pop(host, None)
        try:
            http(url, data=data, use_cache=False, timeout=15)
            out.append((name, True, "OK"))
        except Exception as e:
            out.append((name, False, str(e).split(": ", 1)[-1]))
    try:
        import certifi  # noqa: F401
        out.append(("Certificates", True, "certifi installed"))
    except ImportError:
        out.append(("Certificates", True, "using system certificate store"))
    return out
