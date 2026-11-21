#!/usr/bin/env python3
"""Usage: compare_results.py OUTDIR_A OUTDIR_B  - genes.results expected_count must be identical."""
import csv
import glob
import os
import sys


def load(path):
    with open(path) as fh:
        return {r["gene_id"]: float(r["expected_count"]) for r in csv.DictReader(fh, delimiter="\t")}


def main(argv):
    if len(argv) != 2:
        print("Usage: compare_results.py OUTDIR_A OUTDIR_B")
        return 2
    files = sorted(glob.glob(os.path.join(argv[0], "rsem", "*", "*.genes.results")))
    if not files:
        print("FAIL: no genes.results in " + argv[0])
        return 1
    for fa in files:
        fb = os.path.join(argv[1], "rsem", os.path.basename(os.path.dirname(fa)), os.path.basename(fa))
        a, b = load(fa), load(fb)
        if a.keys() != b.keys() or any(abs(a[g] - b[g]) > 1e-6 for g in a):
            print("FAIL: %s differs between runs" % os.path.basename(fa))
            return 1
    print("ok  %d samples identical" % len(files))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
