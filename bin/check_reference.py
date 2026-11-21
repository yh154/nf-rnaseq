#!/usr/bin/env python3
"""Check a RSEM+STAR reference directory against the GTF, the STAR version and the reads."""
import glob
import os
import re
import sys

VERSION = "0.1.0"
HELP = """\
Usage: check_reference.py --ref_dir DIR --annotation FILE --gtf_contigs FILE --star_version STR
                          [--read_length N] --out_prefix FILE

Check a RSEM+STAR reference (pre-built or just built) before any sample is aligned.

Errors (exit 1):
  * required files missing: <prefix>.grp .ti .seq .transcripts.fa and the STAR index
    (Genome SA SAindex chrName.txt genomeParameters.txt); not exactly one <prefix>.grp
  * STAR index version (versionGenome in genomeParameters.txt) incompatible with the
    STAR of the container (STAR >= 2.7.4a needs an index with versionGenome 2.7.4a)
  * GTF contigs overlap the reference chromosomes by < 90 %
  * a gene of the RSEM reference is absent from the GTF annotation
Warnings (exit 0): sjdbOverhang < read length - 1; GTF genes absent from the reference.

Arguments:
  --ref_dir DIR        reference directory
  --annotation FILE    annotation TSV from build_annotation.py (gene_id column)
  --gtf_contigs FILE   contig names of the GTF exons (build_annotation.py --contigs_out)
  --star_version STR   output of `STAR --version`, e.g. 2.7.10a
  --read_length N      max read length (optional; enables the sjdbOverhang warning)
  --out_prefix FILE    write the RSEM reference prefix (e.g. "ref") here
  -h, --help / --version

Exit codes: 0 ok, 1 reference unusable, 2 usage error.
"""
REQ_SUFFIX = [".grp", ".ti", ".seq", ".transcripts.fa"]
REQ_STAR = ["Genome", "SA", "SAindex", "chrName.txt", "genomeParameters.txt"]


def ver_tuple(v):
    m = re.match(r"(\d+)\.(\d+)\.(\d+)", v)
    return tuple(int(x) for x in m.groups()) if m else None


def read_params(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2:
                d[parts[0]] = parts[1].strip()
    return d


def ti_genes(path):
    """Gene ids from a RSEM .ti file: header 'M type', then 6 lines per transcript."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    try:
        m = int(lines[0].split()[0])
    except (IndexError, ValueError):
        return None
    if len(lines) < 1 + 6 * m:
        return None
    return {lines[1 + 6 * i + 1].split()[0] for i in range(m)}


def check(o):
    errors, warnings = [], []
    d = o["ref_dir"]
    if not os.path.isdir(d):
        return ["reference directory not found: " + d], [], None
    grps = glob.glob(os.path.join(d, "*.grp"))
    if len(grps) != 1:
        return ["expected exactly one <prefix>.grp in %s, found %d" % (d, len(grps))], [], None
    prefix = os.path.basename(grps[0])[:-4]
    for s in REQ_SUFFIX:
        if not os.path.isfile(os.path.join(d, prefix + s)):
            errors.append("missing reference file: %s%s" % (prefix, s))
    for f in REQ_STAR:
        if not os.path.isfile(os.path.join(d, f)):
            errors.append("missing STAR index file: " + f)
    if errors:
        return errors, warnings, prefix

    gp = read_params(os.path.join(d, "genomeParameters.txt"))
    star = ver_tuple(o["star_version"])
    vg = gp.get("versionGenome", "")
    if star is None:
        warnings.append("cannot parse STAR version '%s'; version check skipped" % o["star_version"])
    elif star >= (2, 7, 4) and ver_tuple(vg) != (2, 7, 4):
        errors.append("STAR index versionGenome %s is incompatible with STAR %s (needs 2.7.4a); "
                      "rebuild the reference" % (vg or "?", o["star_version"]))

    if o.get("read_length"):
        try:
            sj = int(gp.get("sjdbOverhang", "0"))
            if sj < o["read_length"] - 1:
                warnings.append("index sjdbOverhang %d < read length - 1 = %d: junction "
                                "sensitivity may drop" % (sj, o["read_length"] - 1))
        except ValueError:
            pass

    with open(os.path.join(d, "chrName.txt")) as fh:
        ref_chr = {l.strip() for l in fh if l.strip()}
    with open(o["gtf_contigs"]) as fh:
        gtf_chr = {l.strip() for l in fh if l.strip()}
    if gtf_chr:
        frac = len(gtf_chr & ref_chr) / len(gtf_chr)
        if frac < 0.9:
            errors.append("only %.0f%% of GTF contigs are in the reference (e.g. GTF %s vs "
                          "reference %s): chr naming mismatch?"
                          % (100 * frac, sorted(gtf_chr)[:2], sorted(ref_chr)[:2]))

    ref_genes = ti_genes(os.path.join(d, prefix + ".ti"))
    with open(o["annotation"]) as fh:
        next(fh, None)
        ann = {l.split("\t")[0] for l in fh if l.strip()}
    if ref_genes is None:
        warnings.append("unrecognised .ti layout: gene-set comparison skipped")
    else:
        missing = ref_genes - ann
        if missing:
            errors.append("%d gene(s) of the reference are not in the GTF (e.g. %s): the GTF "
                          "does not belong to this reference" % (len(missing), sorted(missing)[:3]))
        extra = ann - ref_genes
        if extra:
            warnings.append("%d GTF gene(s) are not in the reference (not quantified)" % len(extra))
    return errors, warnings, prefix


def main(argv):
    o, i = {}, 0
    valued = ("ref_dir", "annotation", "gtf_contigs", "star_version", "read_length", "out_prefix")
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(HELP)
            return 0
        if a == "--version":
            print("check_reference.py " + VERSION)
            return 0
        if a.startswith("--") and a[2:] in valued and i + 1 < len(argv):
            o[a[2:]] = argv[i + 1]
            i += 2
            continue
        sys.stderr.write("ERROR: bad argument: %s\n\n%s" % (a, HELP.split("\n\n")[0] + "\n"))
        return 2
    for r in ("ref_dir", "annotation", "gtf_contigs", "star_version", "out_prefix"):
        if r not in o:
            sys.stderr.write("ERROR: --%s is required\n\n%s" % (r, HELP.split("\n\n")[0] + "\n"))
            return 2
    if "read_length" in o:
        if not o["read_length"].isdigit():
            sys.stderr.write("ERROR: --read_length must be an integer\n")
            return 2
        o["read_length"] = int(o["read_length"])
    errors, warnings, prefix = check(o)
    for w in warnings:
        sys.stderr.write("WARNING: %s\n" % w)
    if errors:
        for e in errors:
            sys.stderr.write("ERROR: %s\n" % e)
        return 1
    with open(o["out_prefix"], "w") as fh:
        fh.write(prefix + "\n")
    print("reference OK: prefix=%s" % prefix)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
