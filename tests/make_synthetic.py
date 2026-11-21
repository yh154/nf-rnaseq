#!/usr/bin/env python3
"""Create the deterministic synthetic test dataset (standard library only)."""
import gzip
import math
import os
import random
import sys

VERSION = "0.1.0"
HELP = """\
Usage: make_synthetic.py [--outdir DIR] [--seed N]

Create a tiny deterministic RNA-seq dataset with known truth:
  genome.fa, genes.gtf (gene_name + gene_biotype, 4 rRNA genes), test.gmt,
  samples_pe.csv / samples_se.csv (2 conditions x 3 replicates, one sample split in
  2 lanes), contrasts.csv, contrast_pairs.txt (same contrast as a plain pair list), truth.tsv (true log2FC = treated/control), fastq/*.fq.gz

Options:
  --outdir DIR   output directory (default: tests/data next to this script)
  --seed N       random seed (default 42)
  -h, --help     show this help
  --version      print version
"""
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def poisson(rng, lam):
    if lam <= 0:
        return 0
    if lam > 30:
        return max(0, int(round(rng.gauss(lam, math.sqrt(lam)))))
    limit, k, p = math.exp(-lam), 0, 1.0
    while True:
        p *= rng.random()
        if p <= limit:
            return k
        k += 1


def build_genes(rng):
    genes, pos, contig = [], 1000, "chr1"
    n_genes, n_rrna = 40, 4
    kinds = ["rRNA" if i in (7, 17, 27, 36) else "protein_coding" for i in range(n_genes)]
    contigs = {"chr1": 0, "chr2": 0}
    for i, kind in enumerate(kinds):
        if i == 24:
            contig, pos = "chr2", 1000
        gid = "G%03d" % (i + 1)
        if kind == "rRNA":
            exons_len = [1800]
        else:
            exons_len = [rng.randint(150, 400) for _ in range(rng.randint(2, 4))]
        exons, p = [], pos
        for k, el in enumerate(exons_len):
            exons.append((p, p + el - 1))
            p += el + (rng.randint(200, 800) if k < len(exons_len) - 1 else 0)
        genes.append({"id": gid, "name": "Gene%d" % (i + 1), "biotype": kind,
                      "contig": contig, "start": exons[0][0], "end": exons[-1][1],
                      "strand": rng.choice("+-"), "exons": exons})
        pos = p + 2000
        contigs[contig] = pos
    assert n_rrna == sum(g["biotype"] == "rRNA" for g in genes)
    return genes, contigs


def write_gtf(path, genes):
    with open(path, "w") as fh:
        for g in genes:
            tid = g["id"] + ".1"
            attr = 'gene_id "%s"; gene_name "%s"; gene_biotype "%s";' % (
                g["id"], g["name"], g["biotype"])
            tattr = 'gene_id "%s"; transcript_id "%s"; gene_name "%s"; gene_biotype "%s";' % (
                g["id"], tid, g["name"], g["biotype"])
            fh.write("\t".join([g["contig"], "synthetic", "gene", str(g["start"]), str(g["end"]),
                                ".", g["strand"], ".", attr]) + "\n")
            fh.write("\t".join([g["contig"], "synthetic", "transcript", str(g["start"]),
                                str(g["end"]), ".", g["strand"], ".", tattr]) + "\n")
            for n, (s, e) in enumerate(g["exons"], 1):
                fh.write("\t".join([g["contig"], "synthetic", "exon", str(s), str(e), ".",
                                    g["strand"], ".", tattr + ' exon_number "%d";' % n]) + "\n")


def transcript_seq(genome, g):
    seq = "".join(genome[g["contig"]][s - 1:e] for s, e in g["exons"])
    return seq if g["strand"] == "+" else revcomp(seq)


def mutate(rng, s, rate=0.005):
    if rng.random() > rate * len(s):
        return s
    i = rng.randrange(len(s))
    return s[:i] + rng.choice([b for b in "ACGT" if b != s[i]]) + s[i + 1:]


def write_fastq(path, records):
    with open(path, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
            for name, seq in records:
                gz.write(("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq))).encode())


def main(argv):
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    seed = 42
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(HELP)
            return 0
        if a == "--version":
            print("make_synthetic.py " + VERSION)
            return 0
        if a in ("--outdir", "--seed") and i + 1 < len(argv):
            if a == "--outdir":
                outdir = argv[i + 1]
            else:
                seed = int(argv[i + 1])
            i += 2
            continue
        sys.stderr.write("ERROR: bad argument %s\n\n%s" % (a, HELP))
        return 2
    rng = random.Random(seed)
    os.makedirs(os.path.join(outdir, "fastq"), exist_ok=True)

    genes, contigs = build_genes(rng)
    genome = {c: "".join(rng.choices("ACGT", k=n + 1000)) for c, n in contigs.items()}
    with open(os.path.join(outdir, "genome.fa"), "w") as fh:
        for c, s in genome.items():
            fh.write(">%s\n" % c)
            for k in range(0, len(s), 60):
                fh.write(s[k:k + 60] + "\n")
    write_gtf(os.path.join(outdir, "genes.gtf"), genes)

    # ---- truth: 8 genes up, 6 down in treated vs control (log2FC = +/-2)
    coding = [g for g in genes if g["biotype"] == "protein_coding"]
    de = rng.sample(coding, 14)
    up, down = de[:8], de[8:]
    base = {}
    for g in genes:
        base[g["id"]] = 1500.0 if g["biotype"] == "rRNA" else math.exp(rng.uniform(math.log(30), math.log(1500)))
    lfc = {g["id"]: 0.0 for g in genes}
    for g in up:
        lfc[g["id"]] = 2.0
    for g in down:
        lfc[g["id"]] = -2.0
    with open(os.path.join(outdir, "truth.tsv"), "w") as fh:
        fh.write("gene_id\tgene_name\tbiotype\tbase_mean\ttrue_log2fc_treated_vs_control\n")
        for g in genes:
            fh.write("%s\t%s\t%s\t%.2f\t%.1f\n" % (g["id"], g["name"], g["biotype"],
                                                    base[g["id"]], lfc[g["id"]]))

    # ---- GMT (gene_name ids): one set holds exactly the up genes, one the down genes
    with open(os.path.join(outdir, "test.gmt"), "w") as fh:
        fh.write("SPIKED_UP\tsynthetic\t" + "\t".join(g["name"] for g in up) + "\n")
        fh.write("SPIKED_DOWN\tsynthetic\t" + "\t".join(g["name"] for g in down) + "\n")
        others = [g for g in coding if g not in de]
        for k in range(1, 7):
            fh.write("RANDOM_SET_%d\tsynthetic\t%s\n" % (
                k, "\t".join(g["name"] for g in rng.sample(others, 8))))

    # ---- reads
    tx = {g["id"]: transcript_seq(genome, g) for g in genes}
    samples = [("control_rep%d" % r, "control", "b%d" % (1 if r != 2 else 2)) for r in (1, 2, 3)]
    samples += [("treated_rep%d" % r, "treated", "b%d" % (1 if r != 2 else 2)) for r in (1, 2, 3)]
    sf = {s[0]: rng.uniform(0.8, 1.2) for s in samples}
    layouts = {"pe": 75, "se": 50}
    fq_dir = os.path.join(outdir, "fastq")
    rows = {"pe": [], "se": []}
    for name, cond, batch in samples:
        frags = []
        for g in genes:
            mu = base[g["id"]] * sf[name] * (2 ** lfc[g["id"]] if cond == "treated" else 1.0)
            phi = 0.05
            lam = rng.gammavariate(1 / phi, mu * phi)
            frags.extend([g["id"]] * poisson(rng, lam))
        rng.shuffle(frags)
        for layout, rl in layouts.items():
            r1, r2 = [], []
            for n, gid in enumerate(frags):
                t = tx[gid]
                flen = min(len(t), max(rl, int(rng.gauss(250, 30))))
                st = rng.randint(0, len(t) - flen)
                f = t[st:st + flen]
                rid = "%s_%s_%06d" % (name, layout, n)
                a = mutate(rng, revcomp(f[-rl:]))     # reverse-stranded: R1 antisense
                r1.append((rid, a))
                if layout == "pe":
                    r2.append((rid, mutate(rng, f[:rl])))
            lanes = [(None, r1, r2)]
            if name == "control_rep3":                 # two lanes for one sample
                h = len(r1) // 2
                lanes = [("L001", r1[:h], r2[:h]), ("L002", r1[h:], r2[h:])]
            for lane, a, b in lanes:
                suffix = "_" + lane if lane else ""
                f1 = "fastq/%s%s_%s_R1.fq.gz" % (name, suffix, layout)
                write_fastq(os.path.join(outdir, f1), a)
                f2 = ""
                if layout == "pe":
                    f2 = "fastq/%s%s_%s_R2.fq.gz" % (name, suffix, layout)
                    write_fastq(os.path.join(outdir, f2), b)
                rows[layout].append([name, f1, f2, "reverse", cond, batch])
    for layout, rr in rows.items():
        with open(os.path.join(outdir, "samples_%s.csv" % layout), "w") as fh:
            fh.write("sample,fastq_1,fastq_2,strandedness,condition,batch\n")
            for r in rr:
                fh.write(",".join(r) + "\n")
    with open(os.path.join(outdir, "contrasts.csv"), "w") as fh:
        fh.write("id,variable,treatment,control\ntreated_vs_control,condition,treated,control\n")
    with open(os.path.join(outdir, "contrast_pairs.txt"), "w") as fh:
        fh.write("# desired contrast pairs: treatment,control\ntreated,control\n")
    print("synthetic data written to %s" % outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
