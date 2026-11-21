#!/usr/bin/env python3
"""Warn when few reads end up in the gene counts (standard library only)."""
import glob
import os
import sys

VERSION = "0.1.0"
HELP = """\
Usage: check_assigned.py --rsem_dir DIR [--min_fraction F] [--out FILE]

Compare, per sample, the sum of the expected counts in <sample>.genes.results with the number of
input reads in the STAR log <sample>.log ("Number of input reads"). A fraction below --min_fraction
is reported as a warning (one line per sample, written to --out and stderr). A low fraction usually
means a wrong --strandedness declaration (which does not always make RSEM fail), the wrong
reference or annotation, or contamination.

Arguments:
  --rsem_dir DIR        directory with <sample>.genes.results and <sample>.log (flat)
  --min_fraction F      warn below this fraction of input reads (default 0.3)
  --out FILE            write the warnings here (empty file when there are none)
  -h, --help / --version

Exit codes: 0 ok (warnings do not change it), 1 unreadable inputs, 2 usage error.

Example:
  check_assigned.py --rsem_dir rsem --min_fraction 0.3 --out assigned_warnings.txt
"""


def usage_error(msg):
    sys.stderr.write("ERROR: %s\n\n%s" % (msg, HELP.split("\n\n")[0] + "\n"))
    return 2


def input_reads(path):
    with open(path) as fh:
        for line in fh:
            if "Number of input reads" in line:
                return int(line.split("|")[1].strip())
    raise ValueError("no 'Number of input reads' line in %s" % path)


def assigned(path):
    total = 0.0
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = header.index("expected_count")
        for line in fh:
            total += float(line.rstrip("\n").split("\t")[col])
    return total


def main(argv):
    o = {"min_fraction": "0.3"}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(HELP)
            return 0
        if a == "--version":
            print("check_assigned.py " + VERSION)
            return 0
        if a in ("--rsem_dir", "--min_fraction", "--out") and i + 1 < len(argv):
            o[a[2:]] = argv[i + 1]
            i += 2
            continue
        return usage_error("bad argument: " + a)
    if "rsem_dir" not in o:
        return usage_error("--rsem_dir is required")
    try:
        minf = float(o["min_fraction"])
    except ValueError:
        return usage_error("--min_fraction must be a number")
    if not 0 <= minf <= 1:
        return usage_error("--min_fraction must be between 0 and 1")
    warnings = []
    try:
        for g in sorted(glob.glob(os.path.join(o["rsem_dir"], "*.genes.results"))):
            sample = os.path.basename(g)[:-len(".genes.results")]
            total = input_reads(os.path.join(o["rsem_dir"], sample + ".log"))
            frac = assigned(g) / total if total else 0.0
            if frac < minf:
                warnings.append("%s: only %.1f%% of the %d input reads are in the gene counts "
                                "(expected > %.0f%%); check --strandedness, the reference and the annotation"
                                % (sample, 100 * frac, total, 100 * minf))
    except (OSError, ValueError, IndexError) as exc:
        sys.stderr.write("ERROR: %s\n" % exc)
        return 1
    for w in warnings:
        sys.stderr.write("WARNING: %s\n" % w)
    if "out" in o:
        with open(o["out"], "w") as fh:
            fh.write("".join(w + "\n" for w in warnings))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
