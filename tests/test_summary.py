#!/usr/bin/env python3
"""Unit tests for bin/summaryHTML.sh (fake FastQC zips and STAR logs, no pipeline needed)."""
import os
import re
import shutil
import subprocess
import tempfile
import unittest
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "bin", "summaryHTML.sh")
MODULES = ["Basic Statistics", "Per base sequence quality", "Per tile sequence quality",
           "Per sequence quality scores", "Per base sequence content", "Per sequence GC content",
           "Per base N content", "Sequence Length Distribution", "Sequence Duplication Levels",
           "Overrepresented sequences", "Adapter Content"]
LOG = """\
                          Number of input reads |\t%d
                   Uniquely mapped reads number |\t%d
        Number of reads mapped to multiple loci |\t%d
        Number of reads mapped to too many loci |\t%d
  Number of reads unmapped: too many mismatches |\t%d
            Number of reads unmapped: too short |\t%d
                Number of reads unmapped: other |\t%d
"""


class TestSummary(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="sum mary ")          # a space in the path on purpose
        self.addCleanup(shutil.rmtree, self.tmp, True)
        for d in ("fastqc", "rsem", "out"):
            os.makedirs(os.path.join(self.tmp, d))

    def p(self, *a):
        return os.path.join(self.tmp, *a)

    def zip(self, base, statuses=None, tile=True):
        statuses = statuses or ["PASS"] * 11
        lines = "".join("%s\t%s\t%s.fq.gz\n" % (s, m, base) for s, m in zip(statuses, MODULES)
                        if tile or "tile" not in m)
        with zipfile.ZipFile(self.p("fastqc", base + "_fastqc.zip"), "w") as z:
            z.writestr("%s_fastqc/summary.txt" % base, lines)

    def sample(self, name, counts=(1000, 900, 60, 10, 5, 15, 10), genes=None):
        with open(self.p("rsem", name + ".log"), "w") as fh:
            fh.write(LOG % counts)
        genes = genes or [("G1", 400.4), ("G2", 99.6), ("G3", 20.0)]
        with open(self.p("rsem", name + ".genes.results"), "w") as fh:
            fh.write("gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM\n")
            for g, c in genes:
                fh.write("%s\tT\t1\t1\t%.2f\t1\t1\n" % (g, c))

    def annotation(self, biotypes):
        with open(self.p("ann.tsv"), "w") as fh:
            fh.write("gene_id\tgene_name\tgene_biotype\n")
            for g, b in biotypes.items():
                fh.write("%s\t%s\t%s\n" % (g, g, b))
        return self.p("ann.tsv")

    def run_script(self, *extra):
        return subprocess.run(["bash", SCRIPT, *extra, self.p("fastqc"), self.p("rsem"), self.p("out")],
                              capture_output=True, text=True)

    def html(self):
        with open(self.p("out", "Summary.html")) as fh:
            return fh.read()

    def test_full_report(self):
        self.zip("A_R1", ["PASS", "WARN", "PASS", "FAIL"] + ["PASS"] * 7)   # index 2 = per-tile (skipped)
        self.zip("A_R2")
        self.sample("A")
        self.sample("B", counts=(200, 150, 0, 0, 20, 20, 10), genes=[("G1", 10), ("G2", 20), ("G3", 5)])
        with open(self.p("order.txt"), "w") as fh:
            fh.write("B\nA\n")
        ann = self.annotation({"G1": "rRNA", "G2": "Mt_rRNA", "G3": "protein_coding"})
        p = self.run_script("--annotation", ann, "--sample_order", self.p("order.txt"),
                            "--footer", "Run <1> & co")
        self.assertEqual(p.returncode, 0, p.stderr)
        h = self.html()
        self.assertEqual(len(re.findall(r"<tr><td>A_R[12]</td>", h)), 2)
        self.assertIn('<td class="warn">WARN</td>', h)
        self.assertIn('<td class="fail">FAIL</td>', h)
        self.assertEqual(len(re.findall(r'<td class="(pass|warn|fail)">', h)), 20)   # 10 modules x 2 files
        self.assertLess(h.index("<td>B</td>"), h.index("<td>A</td>"))              # sample order
        self.assertIn("<tr><td>A</td><td>1000</td><td>900(90.00%)</td><td>70(7.00%)</td>"
                      "<td>30(3.00%)</td><td>500(50.00%)</td></tr>", h)             # rRNA 400.4+99.6
        self.assertIn("<td>30(15.00%)</td>", h)                                    # B: rRNA 30/200
        self.assertIn("Run &lt;1&gt; &amp; co", h)
        self.assertIn("fragments (read pairs)", h)
        self.assertNotIn("@", h)

    def test_single_end_note_and_na_rrna(self):
        self.zip("A")
        self.sample("A")
        for args in ([], ["--annotation", self.annotation({"G1": "NA", "G2": "NA", "G3": "NA"})]):
            p = self.run_script("--single-end", *args)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("<td>NA</td>", self.html())
            self.assertIn("Single-end data: counts are reads.", self.html())

    def test_works_without_tile_module(self):
        self.zip("A", ["PASS"] * 11, tile=False)
        self.sample("A")
        self.assertEqual(self.run_script().returncode, 0)

    def test_inconsistent_log_warns(self):
        self.zip("A")
        self.sample("A", counts=(1000, 900, 60, 10, 5, 15, 5))
        p = self.run_script()
        self.assertEqual(p.returncode, 0)
        self.assertIn("WARNING", p.stderr)

    def test_optional_sections(self):
        self.zip("A")
        self.sample("A")
        with open(self.p("info.txt"), "w") as fh:
            fh.write("Pipeline: nf-rnaseq 0.1.0\nDesign: ~ condition\n")
        with open(self.p("w.txt"), "w") as fh:
            fh.write("one warning\n")
        open(self.p("empty.txt"), "w").close()
        self.assertEqual(self.run_script("--info", self.p("info.txt"), "--warnings", self.p("w.txt")).returncode, 0)
        h = self.html()
        self.assertIn("Run info", h)
        self.assertIn("Validation warnings", h)
        self.assertEqual(self.run_script("--warnings", self.p("empty.txt")).returncode, 0)
        self.assertNotIn("Validation warnings", self.html())

    def test_errors(self):
        self.assertEqual(self.run_script().returncode, 1)                 # no zips
        self.zip("A")
        self.assertEqual(self.run_script().returncode, 1)                 # no logs
        self.sample("A")
        os.remove(self.p("rsem", "A.genes.results"))
        self.assertEqual(self.run_script().returncode, 1)                 # no genes.results
        self.sample("A")
        with open(self.p("rsem", "A.log"), "w") as fh:
            fh.write("garbage\n")
        self.assertEqual(self.run_script().returncode, 1)                 # unparsable log
        with open(self.p("fastqc", "A_fastqc.zip"), "w") as fh:
            fh.write("not a zip")
        self.sample("A")
        self.assertEqual(self.run_script().returncode, 1)                 # corrupt zip

    def test_usage(self):
        self.assertEqual(subprocess.run(["bash", SCRIPT], capture_output=True).returncode, 2)
        self.assertEqual(subprocess.run(["bash", SCRIPT, "--nope"], capture_output=True).returncode, 2)


if __name__ == "__main__":
    unittest.main(verbosity=1)
