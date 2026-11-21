#!/usr/bin/env python3
"""Validate GMT gene-set files against the annotated genes (standard library only)."""
import os
import sys

VERSION = "0.1.0"
HELP = """\
Usage: validate_gmt.py --gmt FILE[,FILE...] --annotation FILE --gmt_id_type gene_name|gene_id
                       [--min_set_size N] [--out FILE]

Check GMT files before GSEA / GSVA so that problems surface early.

Errors (exit 1):
  * a file is unreadable, or two files share a base name (the name labels the outputs)
  * a line has fewer than 3 tab-separated fields (set name, description, >= 1 gene)
  * a gene-set name is repeated inside a file
  * no gene set has at least --min_set_size genes among the annotated genes
Warnings (written to --out and stderr): sets with < 50 % of their genes in the annotation, and a
file where fewer than half of the sets reach --min_set_size annotated genes (usually the wrong
species or --gmt_id_type; matching is exact and case-sensitive).

Arguments:
  --gmt LIST            comma-separated GMT files
  --annotation FILE     annotation TSV (gene_id, gene_name, gene_biotype)
  --gmt_id_type TYPE    which annotation column the GMT genes refer to (default gene_name)
  --min_set_size N      minimum number of annotated genes per usable set (default 5)
  --out FILE            write a one-line-per-file report
  -h, --help / --version

Exit codes: 0 ok, 1 invalid GMT, 2 usage error.

Example:
  validate_gmt.py --gmt h.all.gmt,c2.cp.gmt --annotation annotation.tsv --gmt_id_type gene_name
"""


def usage_error(msg):
    sys.stderr.write("ERROR: %s\n\n%s" % (msg, HELP.split("\n\n")[0] + "\n"))
    return 2


def main(argv):
    o = {"gmt_id_type": "gene_name", "min_set_size": "5"}
    i = 0
    valued = ("gmt", "annotation", "gmt_id_type", "min_set_size", "out")
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(HELP)
            return 0
        if a == "--version":
            print("validate_gmt.py " + VERSION)
            return 0
        if a.startswith("--") and a[2:] in valued and i + 1 < len(argv):
            o[a[2:]] = argv[i + 1]
            i += 2
            continue
        return usage_error("bad argument: " + a)
    for r in ("gmt", "annotation"):
        if r not in o:
            return usage_error("--%s is required" % r)
    if o["gmt_id_type"] not in ("gene_name", "gene_id"):
        return usage_error("--gmt_id_type must be gene_name or gene_id")
    if not o["min_set_size"].isdigit() or int(o["min_set_size"]) < 1:
        return usage_error("--min_set_size must be a positive integer")
    min_size = int(o["min_set_size"])

    known = set()
    with open(o["annotation"]) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = header.index(o["gmt_id_type"])
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > col:
                known.add(f[col])

    errors, warnings, report = [], [], []
    files = [g.strip() for g in o["gmt"].split(",") if g.strip()]
    names = [os.path.basename(f) for f in files]
    if len(set(names)) != len(names):
        errors.append("GMT files must have distinct base names: %s" % names)
    for path in files:
        try:
            with open(path) as fh:
                lines = [l.rstrip("\n") for l in fh if l.strip()]
        except OSError as exc:
            errors.append("cannot read %s: %s" % (path, exc.strerror))
            continue
        seen, usable, nsets = set(), 0, 0
        for n, line in enumerate(lines, 1):
            f = line.split("\t")
            if len(f) < 3 or not f[0] or not any(g.strip() for g in f[2:]):
                errors.append("%s line %d: need set name, description and >= 1 gene (tab-separated)"
                              % (os.path.basename(path), n))
                continue
            nsets += 1
            if f[0] in seen:
                errors.append("%s line %d: repeated gene-set name %s" % (os.path.basename(path), n, f[0]))
            seen.add(f[0])
            genes = {g.strip() for g in f[2:] if g.strip()}
            hit = len(genes & known)
            if hit >= min_size:
                usable += 1
            if hit < 0.5 * len(genes):
                warnings.append("%s: %s has only %d of %d genes in the annotation (%s)"
                                % (os.path.basename(path), f[0], hit, len(genes), o["gmt_id_type"]))
        report.append("%s: %d sets, %d usable (>= %d annotated genes)"
                      % (os.path.basename(path), nsets, usable, min_size))
        if nsets and usable < 0.5 * nsets:
            msg = ("%s: only %d of %d sets have >= %d genes among the annotated %s values "
                   "(wrong species or --gmt_id_type? matching is exact and case-sensitive)"
                   % (os.path.basename(path), usable, nsets, min_size, o["gmt_id_type"]))
            warnings.append(msg)
            report.append("WARNING: " + msg)
        if nsets and usable == 0:
            errors.append("%s: no gene set has >= %d genes among the annotated %s values "
                          "(wrong --gmt_id_type?)" % (os.path.basename(path), min_size, o["gmt_id_type"]))
    for w in warnings:
        sys.stderr.write("WARNING: %s\n" % w)
    if "out" in o:
        with open(o["out"], "w") as fh:
            fh.write("".join(r + "\n" for r in report))
    if errors:
        for e in errors:
            sys.stderr.write("ERROR: %s\n" % e)
        return 1
    print("; ".join(report))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
