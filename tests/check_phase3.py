#!/usr/bin/env python3
"""Assertions on a pipeline output directory after FASTQC + STAR_RSEM (standard library only).

Usage: check_phase3.py OUTDIR SAMPLESHEET [--no-bam] [--samples a,b]
"""
import csv
import glob
import gzip
import math
import os
import struct
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "tests", "data")
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


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


def bam_header(path):
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[-28:] != BGZF_EOF:
        fail("%s: missing BGZF EOF marker (truncated)" % path)
    data = gzip.decompress(raw)
    if data[:4] != b"BAM\1":
        fail("%s: not a BAM" % path)
    (l_text,) = struct.unpack("<i", data[4:8])
    text = data[8:8 + l_text].decode()
    off = 8 + l_text
    (n_ref,) = struct.unpack("<i", data[off:off + 4])
    off += 4
    refs = []
    for _ in range(n_ref):
        (l_name,) = struct.unpack("<i", data[off:off + 4])
        name = data[off + 4:off + 4 + l_name - 1].decode()
        (l_ref,) = struct.unpack("<i", data[off + 4 + l_name:off + 8 + l_name])
        refs.append((name, l_ref))
        off += 8 + l_name
    return text, refs


def fasta_lengths(path):
    out, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = 0
            elif name:
                out[name] += len(line)
    return out


def count_reads(path):
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def read_results(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["gene_id"]: float(r["expected_count"]) for r in rows}


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy)


def main(argv):
    if len(argv) < 2:
        print("Usage: check_phase3.py OUTDIR SAMPLESHEET [--no-bam] [--samples a,b]")
        return 2
    out, sheet = argv[0], argv[1]
    no_bam = "--no-bam" in argv
    only = None
    if "--samples" in argv:
        only = set(argv[argv.index("--samples") + 1].split(","))
    with open(sheet) as fh:
        rows = list(csv.DictReader(fh))
    base = os.path.dirname(os.path.abspath(sheet))
    by_sample = {}
    for r in rows:
        by_sample.setdefault(r["sample"], []).append(r)
    truth, lfc = {}, {}
    with open(os.path.join(DATA, "truth.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            truth[r["gene_id"]] = float(r["base_mean"])
            lfc[r["gene_id"]] = float(r["true_log2fc_treated_vs_control"])
    genome = fasta_lengths(os.path.join(DATA, "genome.fa"))

    for sample, srows in by_sample.items():
        if only and sample not in only:
            continue
        # FastQC: one zip + html per FASTQ file
        for r in srows:
            for col in ("fastq_1", "fastq_2"):
                if r[col]:
                    stem = os.path.basename(r[col])
                    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
                        if stem.endswith(ext):
                            stem = stem[:-len(ext)]
                    for suffix in ("_fastqc.zip", "_fastqc.html"):
                        p = os.path.join(out, "fastqc", stem + suffix)
                        if not os.path.getsize(p) > 0 if os.path.exists(p) else True:
                            fail("missing FastQC output " + p)
        d = rsem_dir(out, sample)
        genes = read_results(os.path.join(d, sample + ".genes.results"))
        if len(genes) != len(truth):
            fail("%s: %d genes in genes.results, expected %d" % (sample, len(genes), len(truth)))
        if not os.path.getsize(os.path.join(d, sample + ".isoforms.results")):
            fail(sample + ": empty isoforms.results")
        if not os.path.isdir(os.path.join(d, sample + ".stat")):
            fail(sample + ": missing RSEM .stat directory")
        # all lanes were used: STAR input reads == FASTQ read (pair) count over lanes
        expected_reads = sum(count_reads(r["fastq_1"] if os.path.isabs(r["fastq_1"])
                                         else os.path.join(base, r["fastq_1"])) for r in srows)
        n_in = None
        with open(os.path.join(d, sample + ".log")) as fh:
            for line in fh:
                if "Number of input reads" in line:
                    n_in = int(line.split("|")[1])
        if n_in != expected_reads:
            fail("%s: STAR input reads %s != FASTQ reads %d (lanes lost or mixed?)"
                 % (sample, n_in, expected_reads))
        # quantification agrees with the simulation truth and (almost) all reads are assigned
        ids = sorted(truth)
        treated = srows[0]["condition"] == "treated"
        expect = {g: truth[g] * (2 ** lfc[g] if treated else 1.0) for g in ids}
        r = pearson([math.log2(genes[g] + 1) for g in ids], [math.log2(expect[g] + 1) for g in ids])
        if r < 0.9:
            fail("%s: log-count correlation with truth %.3f < 0.9 (wrong strandedness?)" % (sample, r))
        assigned = sum(genes.values()) / expected_reads
        if assigned < 0.8:
            fail("%s: only %.0f%% of fragments assigned to genes" % (sample, 100 * assigned))
        # BAM
        bam = os.path.join(out, "bam", sample + ".genome.sorted.bam")
        if no_bam:
            bam_dir = os.path.join(out, "bam")
            if os.path.isdir(bam_dir) and os.listdir(bam_dir):
                fail("--no_bam but bam/ contains files: %s" % os.listdir(bam_dir))
        else:
            text, refs = bam_header(bam)
            if "SO:coordinate" not in text.split("\n")[0]:
                fail(sample + ": BAM is not coordinate-sorted: " + text.split("\n")[0])
            if dict(refs) != genome:
                fail(sample + ": BAM @SQ does not match the FASTA contigs")
            if not os.path.getsize(bam + ".bai") > 0:
                fail(sample + ": empty BAM index")
        for junk in (".STAR.genome.bam", ".transcript.bam", ".genome.bam"):
            for root in (d, os.path.join(out, "bam")):
                if os.path.exists(os.path.join(root, sample + junk)):
                    fail("unexpected intermediate published: " + sample + junk)
        print("ok  %-14s reads=%d assigned=%.0f%% r=%.3f" % (sample, expected_reads, 100 * assigned, r))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
