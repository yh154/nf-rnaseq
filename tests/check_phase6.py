#!/usr/bin/env python3
"""GSEA / GSVA outputs against the simulation truth (standard library only).

Usage: check_phase6.py GSEA_DIR GSVA_DIR [GMT_NAME] [--relaxed-gsva]
--relaxed-gsva: for log2(TPM+1) input, which is compositional (see README): only SPIKED_UP must rank
first and SPIKED_DOWN must be negative with p < 0.05
The synthetic GMT has SPIKED_UP (up in treated), SPIKED_DOWN (down in treated) and 6 random sets.
"""
import csv
import os
import sys


def fail(msg):
    print("FAIL: " + msg)
    sys.exit(1)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    gsea, gsva = argv[0], argv[1]
    relaxed = "--relaxed-gsva" in argv
    argv = [a for a in argv if a != "--relaxed-gsva"]
    name = argv[2] if len(argv) > 2 else "test"
    with open(os.path.join(gsea, "GSEA_treated_vs_control_%s.csv" % name), newline="") as fh:
        g = list(csv.DictReader(fh))
    if [c for c in g[0]] != ["pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge"]:
        fail("GSEA columns: %s" % list(g[0]))
    by = {r["pathway"]: r for r in g}
    up, down = by["SPIKED_UP"], by["SPIKED_DOWN"]
    if not float(up["NES"]) > 0 or not float(down["NES"]) < 0:
        fail("GSEA signs: SPIKED_UP NES %s, SPIKED_DOWN NES %s (positive NES must mean enriched in treated)"
             % (up["NES"], down["NES"]))
    if float(up["padj"]) > 0.05 or float(down["padj"]) > 0.05:
        fail("GSEA: spiked sets not significant (padj %s, %s)" % (up["padj"], down["padj"]))
    nes_abs = sorted(g, key=lambda r: -abs(float(r["NES"])))
    if {nes_abs[0]["pathway"], nes_abs[1]["pathway"]} != {"SPIKED_UP", "SPIKED_DOWN"}:
        fail("GSEA: spiked sets are not the two largest |NES|")
    if [r["pathway"] for r in g][:2] != ["SPIKED_UP", "SPIKED_DOWN"] and \
            {r["pathway"] for r in g[:2]} != {"SPIKED_UP", "SPIKED_DOWN"}:
        fail("GSEA table is not sorted by padj")
    if not up["leadingEdge"] or ";" not in up["leadingEdge"]:
        fail("GSEA leadingEdge not ';'-joined")
    print("ok  GSEA: SPIKED_UP NES %.2f (padj %.1e), SPIKED_DOWN NES %.2f (padj %.1e)"
          % (float(up["NES"]), float(up["padj"]), float(down["NES"]), float(down["padj"])))

    with open(os.path.join(gsva, "GSVA_scores_%s.csv" % name), newline="") as fh:
        sc = list(csv.DictReader(fh))
    if len(sc) != 8:
        fail("GSVA scores: %d gene sets, expected 8" % len(sc))
    with open(os.path.join(gsva, "GSVA_test_treated_vs_control_%s.csv" % name), newline="") as fh:
        t = list(csv.DictReader(fh))
    if list(t[0]) != ["pathway", "mean_treatment", "mean_control", "mean_diff", "pvalue", "padj"]:
        fail("GSVA test columns: %s" % list(t[0]))
    d = {r["pathway"]: r for r in t}
    if relaxed:
        if t[0]["pathway"] != "SPIKED_UP" or not (float(d["SPIKED_DOWN"]["mean_diff"]) < 0
                                                   and float(d["SPIKED_DOWN"]["pvalue"]) < 0.05):
            fail("GSVA (log2 TPM): SPIKED_UP must rank first and SPIKED_DOWN be negative with p < 0.05")
        print("ok  GSVA (log2 TPM, relaxed): SPIKED_UP ranks first, SPIKED_DOWN negative with p < 0.05")
        return 0
    top2 = {t[0]["pathway"], t[1]["pathway"]}
    if top2 != {"SPIKED_UP", "SPIKED_DOWN"}:
        fail("GSVA: top 2 by p-value are %s" % top2)
    if not float(d["SPIKED_UP"]["mean_diff"]) > 0 or not float(d["SPIKED_DOWN"]["mean_diff"]) < 0:
        fail("GSVA mean_diff sign (treatment - control) is wrong")
    for r in t:
        if abs(float(r["mean_diff"]) - (float(r["mean_treatment"]) - float(r["mean_control"]))) > 1e-9:
            fail("GSVA mean_diff != mean_treatment - mean_control for " + r["pathway"])
    print("ok  GSVA: 8 sets scored; top 2 by p-value are the spiked sets with the right signs")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
