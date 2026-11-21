#!/usr/bin/env python3
"""Assertions for the Phase 4 outputs (expression tables + DESeq2), standard library only.

Usage:
  check_phase4.py tables  EXPR_DIR  GENES_RESULTS_DIR  SAMPLESHEET   tables equal per-sample results
  check_phase4.py deseq   DIFFEXP_DIR EXPR_DIR SAMPLESHEET           CSV layout, truth, direction
  check_phase4.py swap    DIFFEXP_A DIFFEXP_B                         swapped contrast flips the sign
  check_phase4.py golden  EXPR_DIR                                    normalized counts vs golden file
"""
import csv
import glob
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "tests", "data")
GOLDEN = os.path.join(ROOT, "tests", "golden", "DESeq2_normalized_counts.txt")
STATS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]


def fail(msg):
    print("FAIL: " + msg)
    sys.exit(1)


def read_tsv(path):
    with open(path, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        return header, list(r)


def samples_of(sheet):
    out = []
    with open(sheet, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["sample"] not in out:
                out.append(r["sample"])
    return out


def num(x):
    return float("nan") if x == "NA" else float(x)


def truth():
    with open(os.path.join(DATA, "truth.tsv")) as fh:
        return {r["gene_id"]: float(r["true_log2fc_treated_vs_control"])
                for r in csv.DictReader(fh, delimiter="\t")}


def cmd_tables(expr, results_dir, sheet):
    samples = samples_of(sheet)
    ref = {}
    for s in samples:
        path = os.path.join(results_dir, s + ".genes.results")
        if not os.path.exists(path):
            path = os.path.join(results_dir, "rsem", s, s + ".genes.results")
        with open(path, newline="") as fh:
            ref[s] = {r["gene_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    orders = []
    for fname, col in (("RSEM_expected_counts.txt", "expected_count"), ("TPM.txt", "TPM"), ("FPKM.txt", "FPKM")):
        header, rows = read_tsv(os.path.join(expr, fname))
        if header != ["gene_id", "gene_biotype", "gene_name"] + samples:
            fail("%s: unexpected header %s" % (fname, header))
        orders.append([r[0] for r in rows])
        for r in rows:
            for s, v in zip(samples, r[3:]):
                if v != ref[s][r[0]][col]:
                    fail("%s: %s/%s = %s but genes.results has %s" % (fname, r[0], s, v, ref[s][r[0]][col]))
    if not (orders[0] == orders[1] == orders[2]):
        fail("gene order differs between the expression tables")
    print("ok  expression tables equal the per-sample genes.results (%d genes)" % len(orders[0]))


def cmd_deseq(diffexp, expr, sheet):
    samples = samples_of(sheet)
    cond = {}
    with open(sheet, newline="") as fh:
        for r in csv.DictReader(fh):
            cond[r["sample"]] = r["condition"]
    t_lfc = truth()
    files = glob.glob(os.path.join(diffexp, "Deseq2_res_*_vs_*.csv"))
    if not files:
        fail("no Deseq2_res_*.csv in " + diffexp)
    ehead, erows = read_tsv(os.path.join(expr, "RSEM_expected_counts.txt"))
    gene_order = [r[0] for r in erows]
    for f in files:
        name = os.path.basename(f)[len("Deseq2_res_"):-len(".csv")]
        trt, ctl = name.split("_vs_")
        with open(f, newline="") as fh:
            rows = list(csv.reader(fh))
        header, body = rows[0], rows[1:]
        sel = [s for s in samples if cond[s] in (trt, ctl)]
        if header != ["gene_id", "gene_name", "gene_biotype"] + sel + STATS:
            fail("%s: unexpected columns %s" % (name, header))
        if [r[0] for r in body] != gene_order:
            fail(name + ": gene order differs from the expression tables")
        d = {r[0]: dict(zip(header, r)) for r in body}
        for g, r in d.items():
            p = num(r["padj"])
            if not (math.isnan(p) or 0 <= p <= 1):
                fail("%s: padj out of range for %s" % (name, g))
        # direction + accuracy against the simulation (treated vs control only)
        if (trt, ctl) == ("treated", "control"):
            spiked = [g for g, v in t_lfc.items() if v != 0]
            nulls = [g for g, v in t_lfc.items() if v == 0]
            hits = [g for g in spiked if num(d[g]["padj"]) < 0.05]
            recall = len(hits) / len(spiked)
            fp = [g for g in nulls if num(d[g]["padj"]) < 0.05]
            if recall < 0.85:
                fail("recall of spiked genes %.2f < 0.85" % recall)
            if len(fp) / len(nulls) > 0.1:
                fail("false positives %d / %d null genes" % (len(fp), len(nulls)))
            wrong_sign = [g for g in spiked if num(d[g]["log2FoldChange"]) * t_lfc[g] <= 0]
            if wrong_sign:
                fail("log2FoldChange has the wrong sign (not treated/control) for %s" % wrong_sign)
            err = sorted(abs(num(d[g]["log2FoldChange"]) - t_lfc[g]) for g in spiked)
            if err[len(err) // 2] > 0.4:
                fail("median |estimated - true log2FC| = %.2f > 0.4" % err[len(err) // 2])
            tr = [s for s in sel if cond[s] == trt]
            cr = [s for s in sel if cond[s] == ctl]
            for g in spiked:
                mt = sum(float(d[g][s]) for s in tr) / len(tr)
                mc = sum(float(d[g][s]) for s in cr) / len(cr)
                if abs(math.log2(mt / mc) - num(d[g]["log2FoldChange"])) > 0.6:
                    fail("%s: log2(mean treated / mean control) disagrees with log2FoldChange" % g)
            print("ok  %s: recall %.2f (%d/%d), %d false positive(s) of %d nulls, median lfc error %.2f"
                  % (name, recall, len(hits), len(spiked), len(fp), len(nulls), err[len(err) // 2]))
        else:
            print("ok  %s layout" % name)


def cmd_swap(a, b):
    fa = glob.glob(os.path.join(a, "Deseq2_res_*.csv"))[0]
    base = os.path.basename(fa)[len("Deseq2_res_"):-4]
    trt, ctl = base.split("_vs_")
    fb = os.path.join(b, "Deseq2_res_%s_vs_%s.csv" % (ctl, trt))
    if not os.path.exists(fb):
        fail("swapping the contrast did not rename the file to " + os.path.basename(fb))
    with open(fa, newline="") as f1, open(fb, newline="") as f2:
        ra, rb = list(csv.DictReader(f1)), list(csv.DictReader(f2))
    for x, y in zip(ra, rb):
        if abs(num(x["log2FoldChange"]) + num(y["log2FoldChange"])) > 1e-9:
            fail("log2FoldChange did not flip sign for " + x["gene_id"])
        px, py = num(x["padj"]), num(y["padj"])
        if not (px == py or (math.isnan(px) and math.isnan(py)) or abs(px - py) < 1e-12):
            fail("padj changed when the contrast was swapped for " + x["gene_id"])
    print("ok  swapping treatment/control renames the file and flips log2FoldChange")


def cmd_golden(expr):
    new_h, new = read_tsv(os.path.join(expr, "DESeq2_normalized_counts.txt"))
    if not os.path.exists(GOLDEN):
        fail("golden file missing: " + GOLDEN)
    old_h, old = read_tsv(GOLDEN)
    if new_h != old_h or [r[0] for r in new] != [r[0] for r in old]:
        fail("normalized counts: header or gene order differs from the golden file")
    for rn, ro in zip(new, old):
        for a, b in zip(rn[3:], ro[3:]):
            if abs(float(a) - float(b)) > 1e-6 * max(1.0, abs(float(b))):
                fail("normalized counts differ from golden at %s: %s vs %s" % (rn[0], a, b))
    print("ok  DESeq2 normalized counts match the golden file (tolerance 1e-6)")


def main(argv):
    cmds = {"tables": (cmd_tables, 3), "deseq": (cmd_deseq, 3), "swap": (cmd_swap, 2), "golden": (cmd_golden, 1)}
    if not argv or argv[0] not in cmds or len(argv) - 1 != cmds[argv[0]][1]:
        print(__doc__)
        return 2
    cmds[argv[0]][0](*argv[1:])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
