#!/usr/bin/env python3
"""Unit tests for bin/validate_gmt.py."""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "bin", "validate_gmt.py")


class TestGmt(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.w("ann.tsv", "gene_id\tgene_name\tgene_biotype\n" +
               "".join("G%d\tName%d\tprotein_coding\n" % (i, i) for i in range(1, 21)))

    def w(self, name, text):
        with open(os.path.join(self.tmp, name), "w") as fh:
            fh.write(text)
        return os.path.join(self.tmp, name)

    def run_v(self, *gmts, extra=()):
        return subprocess.run([sys.executable, "-I", SCRIPT, "--gmt", ",".join(gmts),
                               "--annotation", os.path.join(self.tmp, "ann.tsv"), *extra],
                              capture_output=True, text=True)

    def gmt(self, name, *sets):
        return self.w(name, "".join("%s\tdesc\t%s\n" % (n, "\t".join(g)) for n, g in sets))

    def names(self, a, b):
        return ["Name%d" % i for i in range(a, b + 1)]

    def test_mostly_unmatched_sets_warn(self):
        out = os.path.join(self.tmp, "r.txt")
        p = self.run_v(self.gmt("a.gmt", ("S1", self.names(1, 8)), ("S2", ["X1", "X2", "X3", "X4", "X5"]),
                                ("S3", ["Y1", "Y2", "Y3", "Y4", "Y5"])),
                       extra=("--min_set_size", "5", "--out", out))
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("only 1 of 3 sets", p.stderr)
        self.assertIn("WARNING: a.gmt: only 1 of 3 sets", open(out).read())

    def test_valid(self):
        p = self.run_v(self.gmt("a.gmt", ("S1", self.names(1, 8)), ("S2", self.names(9, 20))),
                       extra=("--min_set_size", "5", "--out", os.path.join(self.tmp, "r.txt")))
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("2 sets, 2 usable", open(os.path.join(self.tmp, "r.txt")).read())

    def test_gene_id_type(self):
        g = self.gmt("a.gmt", ("S1", ["G%d" % i for i in range(1, 9)]))
        self.assertEqual(self.run_v(g, extra=("--min_set_size", "5", "--gmt_id_type", "gene_id")).returncode, 0)
        p = self.run_v(g, extra=("--min_set_size", "5"))          # default gene_name: nothing matches
        self.assertEqual(p.returncode, 1)
        self.assertIn("wrong --gmt_id_type", p.stderr)

    def test_no_usable_set_when_too_small(self):
        p = self.run_v(self.gmt("a.gmt", ("S1", self.names(1, 4))), extra=("--min_set_size", "5"))
        self.assertEqual(p.returncode, 1)

    def test_malformed_line(self):
        p = self.run_v(self.w("a.gmt", "S1\tonly-two-fields\n"))
        self.assertEqual(p.returncode, 1)
        self.assertIn("line 1", p.stderr)

    def test_duplicate_set_name(self):
        p = self.run_v(self.gmt("a.gmt", ("S1", self.names(1, 8)), ("S1", self.names(9, 16))),
                       extra=("--min_set_size", "5"))
        self.assertEqual(p.returncode, 1)
        self.assertIn("repeated gene-set name", p.stderr)

    def test_missing_file_and_duplicate_basenames(self):
        self.assertEqual(self.run_v(os.path.join(self.tmp, "nope.gmt")).returncode, 1)
        os.makedirs(os.path.join(self.tmp, "x"))
        g1 = self.gmt("a.gmt", ("S1", self.names(1, 8)))
        g2 = self.w("x/a.gmt", "S1\td\t" + "\t".join(self.names(1, 8)) + "\n")
        p = self.run_v(g1, g2, extra=("--min_set_size", "5"))
        self.assertEqual(p.returncode, 1)
        self.assertIn("distinct base names", p.stderr)

    def test_warning_for_unmatched_genes(self):
        p = self.run_v(self.gmt("a.gmt", ("S1", self.names(1, 8)), ("S2", ["X1", "X2", "X3"] + self.names(1, 1))),
                       extra=("--min_set_size", "5"))
        self.assertEqual(p.returncode, 0)
        self.assertIn("WARNING", p.stderr)

    def test_usage(self):
        self.assertEqual(subprocess.run([sys.executable, "-I", SCRIPT], capture_output=True).returncode, 2)
        self.assertEqual(self.run_v("x.gmt", extra=("--gmt_id_type", "symbol")).returncode, 2)


if __name__ == "__main__":
    unittest.main(verbosity=1)
