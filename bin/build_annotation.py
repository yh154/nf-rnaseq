#!/usr/bin/env python3
"""Build the gene annotation table (gene_id, gene_name, gene_biotype) from a GTF."""
import gzip
import re
import sys

VERSION = "0.1.0"
HELP = """\
Usage: build_annotation.py --gtf FILE --out FILE [--contigs_out FILE]

Build the gene annotation table from a GTF (plain or .gz; standard library only).

Only genes with at least one exon feature are listed (these are the genes RSEM quantifies).
  gene_name     attribute gene_name; falls back to gene_id when absent
  gene_biotype  attribute gene_biotype (Ensembl) or gene_type (GENCODE); NA when absent

Options:
  --gtf FILE          input GTF
  --out FILE          output TSV: gene_id, gene_name, gene_biotype (first-seen gene order)
  --contigs_out FILE  also write the contig names that carry exon features, one per line
  -h, --help          show this help and exit 0
  --version           print the version and exit 0

Exit codes: 0 ok, 1 unusable GTF, 2 usage error.

Example:
  build_annotation.py --gtf genome.gtf --out annotation.tsv --contigs_out gtf_contigs.txt
"""
ATTR = re.compile(r'([A-Za-z_][\w.]*)\s+"([^"]*)"')


def parse(path):
    opener = gzip.open if path.endswith(".gz") else open
    genes, contigs = {}, []
    seen_contig = set()
    with opener(path, "rt") as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                raise ValueError("line %d: expected 9 tab-separated fields" % n)
            attrs = dict(ATTR.findall(f[8]))
            gid = attrs.get("gene_id")
            if not gid:
                continue
            g = genes.setdefault(gid, {"name": "", "biotype": "", "exon": False})
            if attrs.get("gene_name") and not g["name"]:
                g["name"] = attrs["gene_name"]
            bt = attrs.get("gene_biotype") or attrs.get("gene_type")
            if bt and not g["biotype"]:
                g["biotype"] = bt
            if f[2] == "exon":
                g["exon"] = True
                if f[0] not in seen_contig:
                    seen_contig.add(f[0])
                    contigs.append(f[0])
    return genes, contigs


def main(argv):
    opts, i = {}, 0
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(HELP)
            return 0
        if a == "--version":
            print("build_annotation.py " + VERSION)
            return 0
        if a in ("--gtf", "--out", "--contigs_out") and i + 1 < len(argv):
            opts[a[2:]] = argv[i + 1]
            i += 2
            continue
        sys.stderr.write("ERROR: bad argument: %s\n\n%s" % (a, HELP.split("\n\n")[0] + "\n"))
        return 2
    if "gtf" not in opts or "out" not in opts:
        sys.stderr.write("ERROR: --gtf and --out are required\n\n%s" % (HELP.split("\n\n")[0] + "\n"))
        return 2
    try:
        genes, contigs = parse(opts["gtf"])
    except (OSError, ValueError) as exc:
        sys.stderr.write("ERROR: cannot read GTF %s: %s\n" % (opts["gtf"], exc))
        return 1
    rows = [(gid, g) for gid, g in genes.items() if g["exon"]]
    if not rows:
        sys.stderr.write("ERROR: no gene with exon features and gene_id in %s\n" % opts["gtf"])
        return 1
    with open(opts["out"], "w") as fh:
        fh.write("gene_id\tgene_name\tgene_biotype\n")
        for gid, g in rows:
            fh.write("%s\t%s\t%s\n" % (gid, g["name"] or gid, g["biotype"] or "NA"))
    if "contigs_out" in opts:
        with open(opts["contigs_out"], "w") as fh:
            fh.write("".join(c + "\n" for c in contigs))
    nbt = sum(1 for _, g in rows if g["biotype"])
    print("annotation: %d genes, %d with biotype, %d contigs" % (len(rows), nbt, len(contigs)))
    if nbt == 0:
        sys.stderr.write("WARNING: no gene_biotype/gene_type attribute found; biotype is NA "
                         "(Summary rRNA column will show NA)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
