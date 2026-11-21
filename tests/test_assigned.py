#!/usr/bin/env python3
"""Unit tests for bin/check_assigned.py."""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "bin", "check_assigned.py")


class TestAssigned(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def sample(self, name, reads, counts):
        with open(os.path.join(self.tmp, name + ".log"), "w") as fh:
            fh.write("                          Number of input reads |\t%d\n" % reads)
        with open(os.path.join(self.tmp, name + ".genes.results"), "w") as fh:
            fh.write("gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM\n")
            for i, c in enumerate(counts):
                fh.write("g%d\tt%d\t1000\t900\t%.2f\t1\t1\n" % (i, i, c))

    def run_c(self, *extra):
        out = os.path.join(self.tmp, "w.txt")
        p = subprocess.run([sys.executable, "-I", SCRIPT, "--rsem_dir", self.tmp, "--out", out, *extra],
                           capture_output=True, text=True)
        return p, open(out).read() if os.path.exists(out) else None

    def test_good_and_low(self):
        self.sample("good", 1000, [400, 300])
        self.sample("bad", 1000, [50, 40])
        p, w = self.run_c()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("bad: only 9.0%", w)
        self.assertNotIn("good", w)
        self.assertEqual(len(w.splitlines()), 1)

    def test_no_warning_gives_empty_file(self):
        self.sample("good", 1000, [900])
        p, w = self.run_c()
        self.assertEqual((p.returncode, w), (0, ""))

    def test_threshold_option(self):
        self.sample("s", 1000, [400])
        p, w = self.run_c("--min_fraction", "0.5")
        self.assertIn("s: only 40.0%", w)

    def test_missing_log_is_error(self):
        self.sample("s", 1000, [400])
        os.remove(os.path.join(self.tmp, "s.log"))
        p, w = self.run_c()
        self.assertEqual(p.returncode, 1)

    def test_bad_args(self):
        p = subprocess.run([sys.executable, "-I", SCRIPT, "--min_fraction", "2", "--rsem_dir", self.tmp],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 2)


if __name__ == "__main__":
    unittest.main()
