#!/usr/bin/env python3
"""Merge per-sample RSEM genes.results into gene-level expression tables."""
import csv
import os
import sys

VERSION = "0.1.0"
HELP = """\
Usage: merge_rsem.py --samplesheet FILE --annotation FILE --outdir DIR <sample>.genes.results ...

Merge RSEM gene-level results into three tables (standard library only):
  RSEM_expected_counts.txt   expected_count
  TPM.txt                    TPM
  FPKM.txt                   FPKM
Columns: gene_id, gene_biotype, gene_name, <one column per sample>. Gene order follows the
annotation table, sample order follows the samplesheet (first appearance); the order is
identical in all tables. Values are copied verbatim from genes.results.

Each results file must be named <sample>.genes.results and every sample of the samplesheet
must have one. All files must contain exactly the genes of the annotation.

Options:
  --samplesheet FILE   validated samplesheet (sample column)
  --annotation FILE    annotation TSV: gene_id, gene_name, gene_biotype
  --outdir DIR         output directory (created)
  -h, --help           show this help and exit 0
  --version            print the version and exit 0

Exit codes: 0 ok, 1 inconsistent inputs, 2 usage error.

Example:
  merge_rsem.py --samplesheet samplesheet.valid.csv --annotation annotation.tsv \\
      --outdir Expression *.genes.results
"""
TABLES = [("RSEM_expected_counts.txt", "expected_count"), ("TPM.txt", "TPM"), ("FPKM.txt", "FPKM")]


def usage_error(msg):
    sys.stderr.write("ERROR: %s\n\n%s" % (msg, HELP.split("\n\n")[0] + "\n"))
    return 2


def read_results(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["gene_id"]: r for r in rows}


def main(argv):
    o, files, i = {}, [], 0
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(HELP)
            return 0
        if a == "--version":
            print("merge_rsem.py " + VERSION)
            return 0
        if a in ("--samplesheet", "--annotation", "--outdir") and i + 1 < len(argv):
            o[a[2:]] = argv[i + 1]
            i += 2
            continue
        if a.startswith("-"):
            return usage_error("bad argument: " + a)
        files.append(a)
        i += 1
    for r in ("samplesheet", "annotation", "outdir"):
        if r not in o:
            return usage_error("--%s is required" % r)
    if not files:
        return usage_error("no genes.results files given")

    samples = []
    with open(o["samplesheet"], newline="") as fh:
        for r in csv.DictReader(fh):
            if r["sample"] not in samples:
                samples.append(r["sample"])
    by_sample = {}
    for f in files:
        name = os.path.basename(f)
        if not name.endswith(".genes.results"):
            sys.stderr.write("ERROR: not a <sample>.genes.results file: %s\n" % f)
            return 1
        by_sample[name[:-len(".genes.results")]] = f
    missing = [s for s in samples if s not in by_sample]
    extra = [s for s in by_sample if s not in samples]
    if missing or extra:
        sys.stderr.write("ERROR: samples and results files disagree (missing results: %s; "
                         "unexpected files: %s)\n" % (missing, extra))
        return 1

    with open(o["annotation"], newline="") as fh:
        ann = list(csv.DictReader(fh, delimiter="\t"))
    genes = [a["gene_id"] for a in ann]
    res = {s: read_results(by_sample[s]) for s in samples}
    for s in samples:
        if set(res[s]) != set(genes):
            only_res = sorted(set(res[s]) - set(genes))[:3]
            only_ann = sorted(set(genes) - set(res[s]))[:3]
            sys.stderr.write("ERROR: genes of %s.genes.results differ from the annotation "
                             "(only in results: %s; only in annotation: %s)\n"
                             % (s, only_res, only_ann))
            return 1

    os.makedirs(o["outdir"], exist_ok=True)
    for fname, col in TABLES:
        with open(os.path.join(o["outdir"], fname), "w") as fh:
            fh.write("\t".join(["gene_id", "gene_biotype", "gene_name"] + samples) + "\n")
            for a in ann:
                fh.write("\t".join([a["gene_id"], a["gene_biotype"], a["gene_name"]]
                                   + [res[s][a["gene_id"]][col] for s in samples]) + "\n")
    print("merged %d samples x %d genes" % (len(samples), len(genes)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
