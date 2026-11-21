#!/usr/bin/env python3
"""Unit tests for bin/merge_rsem.py."""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "bin", "merge_rsem.py")
HEADER = "gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM\n"


class TestMerge(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.w("ss.csv", "sample,fastq_1,fastq_2,strandedness,condition\nB,x,y,none,c\nA,x,y,none,c\nA,x2,y2,none,c\n")
        self.w("ann.tsv", "gene_id\tgene_name\tgene_biotype\nG2\tTwo\trRNA\nG1\tOne\tprotein_coding\n")

    def w(self, name, text):
        with open(os.path.join(self.tmp, name), "w") as fh:
            fh.write(text)
        return os.path.join(self.tmp, name)

    def results(self, sample, genes=("G1", "G2"), scale=1):
        rows = "".join("%s\t%s.1\t100\t50\t%d.00\t%d.50\t%d.25\n" % (g, g, scale * (i + 1), scale * (i + 1), scale * (i + 1))
                       for i, g in enumerate(genes))
        return self.w(sample + ".genes.results", HEADER + rows)

    def run_merge(self, *files, ss="ss.csv"):
        return subprocess.run([sys.executable, "-I", SCRIPT, "--samplesheet", os.path.join(self.tmp, ss),
                               "--annotation", os.path.join(self.tmp, "ann.tsv"),
                               "--outdir", os.path.join(self.tmp, "out"), *files],
                              capture_output=True, text=True)

    def read(self, name):
        with open(os.path.join(self.tmp, "out", name)) as fh:
            return fh.read().splitlines()

    def test_order_and_values(self):
        p = self.run_merge(self.results("A"), self.results("B", scale=10))
        self.assertEqual(p.returncode, 0, p.stderr)
        # sample order = samplesheet (B, A); gene order = annotation (G2, G1); values verbatim
        self.assertEqual(self.read("RSEM_expected_counts.txt"),
                         ["gene_id\tgene_biotype\tgene_name\tB\tA",
                          "G2\trRNA\tTwo\t20.00\t2.00", "G1\tprotein_coding\tOne\t10.00\t1.00"])
        self.assertEqual(self.read("TPM.txt")[1], "G2\trRNA\tTwo\t20.50\t2.50")
        self.assertEqual(self.read("FPKM.txt")[2], "G1\tprotein_coding\tOne\t10.25\t1.25")

    def test_missing_sample(self):
        p = self.run_merge(self.results("A"))
        self.assertEqual(p.returncode, 1)
        self.assertIn("missing results", p.stderr)

    def test_gene_set_mismatch(self):
        p = self.run_merge(self.results("A"), self.results("B", genes=("G1", "G3")))
        self.assertEqual(p.returncode, 1)
        self.assertIn("differ from the annotation", p.stderr)

    def test_bad_file_name_and_usage(self):
        self.assertEqual(self.run_merge(self.w("x.txt", "a"), self.results("A")).returncode, 1)
        self.assertEqual(subprocess.run([sys.executable, "-I", SCRIPT], capture_output=True).returncode, 2)


if __name__ == "__main__":
    unittest.main(verbosity=1)
