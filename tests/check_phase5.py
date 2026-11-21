#!/usr/bin/env python3
"""Assertions for Summary.html and the plots PDF (standard library only).

Usage: check_phase5.py OUTDIR SAMPLESHEET N_CONTRASTS RUN_DATE [EXTRA_PAGES]
(EXTRA_PAGES = GSEA + GSVA pages, one each per contrast x GMT file)
"""
import csv
import glob
import os
import re
import sys
from html.parser import HTMLParser
from xml.dom import minidom

PCT = re.compile(r"^(\d+)\((\d+\.\d\d)%\)$")


def fail(msg):
    print("FAIL: " + msg)
    sys.exit(1)


def rsem_dir(out, sample):
    """Task work directory of STAR_RSEM (sample): the per-sample RSEM files are not published."""
    work = os.environ.get("NXF_WORK", "work")
    with open(os.path.join(out, "pipeline_info", "trace.txt")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["name"] == "STAR_RSEM (%s)" % sample and r["status"] in ("COMPLETED", "CACHED"):
                hits = glob.glob(os.path.join(work, r["hash"] + "*"))
                if len(hits) == 1:
                    return hits[0]
    fail("no STAR_RSEM work directory found for " + sample)


class Tables(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tables, self.row, self.cell, self.cls, self.h = [], None, None, None, []
        self.in_h = None
        self.text = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "table":
            self.tables.append([])
        elif tag == "tr":
            self.row = []
        elif tag in ("td", "th"):
            self.cell, self.cls = "", a.get("class", "")
        elif tag in ("h1", "h2"):
            self.in_h = ""

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self.row is not None:
            self.row.append((self.cell, self.cls, tag))
            self.cell = None
        elif tag == "tr" and self.row is not None:
            self.tables[-1].append(self.row)
            self.row = None
        elif tag in ("h1", "h2"):
            self.h.append(self.in_h)
            self.in_h = None

    def handle_data(self, data):
        if self.cell is not None:
            self.cell += data
        if self.in_h is not None:
            self.in_h += data
        self.text.append(data)


def main(argv):
    if len(argv) not in (4, 5):
        print(__doc__)
        return 2
    out, sheet, n_contrasts, run_date = argv[0], argv[1], int(argv[2]), argv[3]
    extra = int(argv[4]) if len(argv) == 5 else 0
    path = os.path.join(out, "Summary.html")
    raw = open(path).read()
    minidom.parseString(raw)                       # well-formed
    if "@" in raw:
        fail("Summary.html contains an '@' (personal e-mail in the footer?)")
    p = Tables()
    p.feed(raw)
    if p.h[:3] != ["Summary Report", "Read Quality Summary", "Mapping Statistics"]:
        fail("unexpected headings %s" % p.h)
    if "#000080" not in raw or "monospace" not in raw:
        fail("reference style block missing")
    fq, mp = p.tables[0], p.tables[1]
    with open(sheet, newline="") as fh:
        rows = list(csv.DictReader(fh))
    samples = []
    for r in rows:
        if r["sample"] not in samples:
            samples.append(r["sample"])
    n_fastq = sum(1 for r in rows for c in ("fastq_1", "fastq_2") if r[c])
    if len(fq) - 1 != n_fastq:
        fail("section 1 has %d rows, expected one per FASTQ (%d)" % (len(fq) - 1, n_fastq))
    for r in fq[1:]:
        if len(r) != 11 or any(c[0] not in ("PASS", "WARN", "FAIL") for c in r[1:]):
            fail("bad FastQC row %s" % [c[0] for c in r])
        for text, cls, _ in r[1:]:
            if cls != text.lower():
                fail("cell colour class %r for %s" % (cls, text))
        if r[0][0].endswith("_fastqc") or r[0][0].endswith(".gz"):
            fail("label not cleaned: " + r[0][0])
    if [r[0][0] for r in mp[1:]] != samples:
        fail("section 2 rows %s != samples %s" % ([r[0][0] for r in mp[1:]], samples))
    hdr = [c[0] for c in mp[0]]
    if hdr != ["Library", "Total Reads", "Uniquely Mapped", "Multi-reads", "Unmapped Reads", "rRNA Reads"]:
        fail("mapping header " + str(hdr))
    # annotation: rRNA genes
    ann = {}
    with open(os.path.join(out, "reference", "annotation.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ann[r["gene_id"]] = r["gene_biotype"]
    for r in mp[1:]:
        s = r[0][0]
        if not re.fullmatch(r"\d+", r[1][0]):
            fail("%s: Total Reads not an integer: %s" % (s, r[1][0]))
        total = int(r[1][0])
        vals = []
        for c in r[2:]:
            m = PCT.match(c[0])
            if not m:
                fail("%s: %r is not N(P.PP%%)" % (s, c[0]))
            n, pct = int(m.group(1)), float(m.group(2))
            if abs(pct - 100.0 * n / total) > 0.0051:
                fail("%s: percentage %s inconsistent with %d/%d" % (s, pct, n, total))
            vals.append(n)
        if vals[0] + vals[1] + vals[2] != total:
            fail("%s: unique+multi+unmapped = %d != total %d" % (s, sum(vals[:3]), total))
        with open(os.path.join(rsem_dir(out, s), s + ".genes.results")) as fh:
            rr = sum(float(g["expected_count"]) for g in csv.DictReader(fh, delimiter="\t")
                     if ann.get(g["gene_id"]) in ("rRNA", "Mt_rRNA", "rRNA_pseudogene"))
        if abs(vals[3] - round(rr)) > 1:
            fail("%s: rRNA reads %d != sum of rRNA expected counts %.1f" % (s, vals[3], rr))
    # plots PDF: name and page count (1 global PCA + 4 per contrast)
    pdfs = glob.glob(os.path.join(out, "plot", "plots_*.pdf"))
    if len(pdfs) != 1 or os.path.basename(pdfs[0]) != "plots_%s.pdf" % run_date:
        fail("expected plot/plots_%s.pdf, found %s" % (run_date, [os.path.basename(x) for x in pdfs]))
    data = open(pdfs[0], "rb").read()
    pages = len(re.findall(rb"/Type\s*/Page\b(?!s)", data))
    if pages != 1 + 4 * n_contrasts + extra:
        fail("PDF has %d pages, expected %d" % (pages, 1 + 4 * n_contrasts + extra))
    print("ok  Summary.html (%d FASTQ rows, %d samples) and %s (%d pages)"
          % (len(fq) - 1, len(samples), os.path.basename(pdfs[0]), pages))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
