#!/usr/bin/env python3
"""Unit tests for Phase 2 scripts: build_annotation, max_read_length, check_reference, hash_inputs."""
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
DATA = os.path.join(ROOT, "tests", "data")


def run(script, *args):
    p = subprocess.run([sys.executable, "-I", os.path.join(BIN, script), *args],
                       capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


class Tmp(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def p(self, name):
        return os.path.join(self.tmp, name)

    def write(self, name, text, mode="w"):
        with open(self.p(name), mode) as fh:
            fh.write(text)
        return self.p(name)

    def read(self, name):
        with open(self.p(name)) as fh:
            return fh.read()


def exon(gene_attrs, contig="1", gene="G1"):
    return "%s\tsrc\texon\t1\t100\t.\t+\t.\tgene_id \"%s\"; transcript_id \"T\"; %s\n" % (
        contig, gene, gene_attrs)


class TestAnnotation(Tmp):
    def build(self, gtf_text, name="a.gtf"):
        g = self.write(name, gtf_text, "wt") if not name.endswith(".gz") else None
        if g is None:
            g = self.p(name)
            with gzip.open(g, "wt") as fh:
                fh.write(gtf_text)
        return run("build_annotation.py", "--gtf", g, "--out", self.p("a.tsv"),
                   "--contigs_out", self.p("c.txt"))

    def test_ensembl_biotype_and_name(self):
        rc, out, _ = self.build(exon('gene_name "Actb"; gene_biotype "protein_coding";'))
        self.assertEqual(rc, 0)
        self.assertEqual(self.read("a.tsv").splitlines()[1], "G1\tActb\tprotein_coding")
        self.assertEqual(self.read("c.txt"), "1\n")

    def test_gencode_gene_type(self):
        self.build(exon('gene_name "X"; gene_type "rRNA";'))
        self.assertTrue(self.read("a.tsv").splitlines()[1].endswith("\trRNA"))

    def test_missing_name_falls_back_and_biotype_na(self):
        rc, _, err = self.build(exon(""))
        self.assertEqual(rc, 0)
        self.assertEqual(self.read("a.tsv").splitlines()[1], "G1\tG1\tNA")
        self.assertIn("WARNING", err)

    def test_genes_without_exon_skipped_and_order_kept(self):
        txt = ('1\ts\tgene\t1\t9\t.\t+\t.\tgene_id "G0"; gene_name "Zero";\n'
               + exon('gene_name "B";', gene="G2") + exon('gene_name "A";', gene="G1"))
        self.build(txt)
        self.assertEqual([l.split("\t")[0] for l in self.read("a.tsv").splitlines()[1:]],
                         ["G2", "G1"])

    def test_gzip_gtf(self):
        rc, _, _ = self.build(exon('gene_name "A";'), "a.gtf.gz")
        self.assertEqual(rc, 0)

    def test_name_taken_from_any_feature(self):
        txt = ('1\ts\tgene\t1\t9\t.\t+\t.\tgene_id "G1"; gene_name "Late"; gene_biotype "lncRNA";\n'
               + exon("", gene="G1"))
        self.build(txt)
        self.assertEqual(self.read("a.tsv").splitlines()[1], "G1\tLate\tlncRNA")

    def test_bad_gtf_exit_1(self):
        self.assertEqual(self.build("not a gtf\n")[0], 1)
        self.assertEqual(self.build("# only comments\n")[0], 1)

    def test_usage_error(self):
        self.assertEqual(run("build_annotation.py")[0], 2)

    def test_synthetic_gtf(self):
        rc, out, _ = run("build_annotation.py", "--gtf", os.path.join(DATA, "genes.gtf"),
                         "--out", self.p("a.tsv"))
        self.assertEqual(rc, 0)
        rows = [l.split("\t") for l in self.read("a.tsv").splitlines()[1:]]
        self.assertEqual(len(rows), 40)
        self.assertEqual(sum(1 for r in rows if r[2] == "rRNA"), 4)


class TestReadLength(Tmp):
    def fq(self, name, lengths, gz=True):
        text = "".join("@r%d\n%s\n+\n%s\n" % (i, "A" * n, "I" * n) for i, n in enumerate(lengths))
        path = self.p(name)
        if gz:
            with gzip.open(path, "wt") as fh:
                fh.write(text)
        else:
            with open(path, "w") as fh:
                fh.write(text)
        return path

    def test_max_over_mates(self):
        self.assertEqual(run("max_read_length.py", self.fq("a.fq.gz", [50, 60]),
                             self.fq("b.fq.gz", [75, 30]))[1].strip(), "75")

    def test_plain_and_n_reads(self):
        f = self.fq("a.fq", [50, 51, 99], gz=False)
        self.assertEqual(run("max_read_length.py", "--n_reads", "2", f)[1].strip(), "51")

    def test_synthetic(self):
        self.assertEqual(run("max_read_length.py", os.path.join(
            DATA, "fastq", "control_rep1_pe_R1.fq.gz"))[1].strip(), "75")
        self.assertEqual(run("max_read_length.py", os.path.join(
            DATA, "fastq", "control_rep1_se_R1.fq.gz"))[1].strip(), "50")

    def test_errors(self):
        self.assertEqual(run("max_read_length.py")[0], 2)
        self.assertEqual(run("max_read_length.py", self.p("missing.fq"))[0], 1)
        self.assertEqual(run("max_read_length.py", self.write("e.fq", ""))[0], 1)


class TestCheckReference(Tmp):
    def make_ref(self, genes=("G1", "G2"), version="2.7.4a", overhang=74, chrs=("1", "2")):
        d = self.p("ref")
        shutil.rmtree(d, True)
        os.makedirs(d)
        for f in ("ref.grp", "ref.seq", "ref.transcripts.fa", "Genome", "SA", "SAindex"):
            open(os.path.join(d, f), "w").close()
        with open(os.path.join(d, "ref.ti"), "w") as fh:
            fh.write("%d 0\n" % len(genes))
            for g in genes:
                fh.write("%s.1\n%s\t%s\n0 1 2\n+\n10 0\n1\n" % (g, g, g))
        with open(os.path.join(d, "chrName.txt"), "w") as fh:
            fh.write("".join(c + "\n" for c in chrs))
        with open(os.path.join(d, "genomeParameters.txt"), "w") as fh:
            fh.write("### STAR\nversionGenome\t%s\nsjdbOverhang\t%d\n" % (version, overhang))
        return d

    def check(self, d, ann=("G1", "G2"), contigs=("1", "2"), star="2.7.10a", rl=None):
        self.write("ann.tsv", "gene_id\tgene_name\tgene_biotype\n" + "".join("%s\t%s\tNA\n" % (g, g) for g in ann))
        self.write("ctg.txt", "".join(c + "\n" for c in contigs))
        args = ["--ref_dir", d, "--annotation", self.p("ann.tsv"), "--gtf_contigs", self.p("ctg.txt"),
                "--star_version", star, "--out_prefix", self.p("prefix.txt")]
        if rl:
            args += ["--read_length", str(rl)]
        return run("check_reference.py", *args)

    def test_ok_and_prefix(self):
        rc, out, err = self.check(self.make_ref())
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.read("prefix.txt").strip(), "ref")

    def test_missing_file(self):
        d = self.make_ref()
        os.remove(os.path.join(d, "SAindex"))
        self.assertIn("missing STAR index file", self.check(d)[2])

    def test_missing_rsem_file(self):
        d = self.make_ref()
        os.remove(os.path.join(d, "ref.seq"))
        rc, _, err = self.check(d)
        self.assertEqual(rc, 1)
        self.assertIn("ref.seq", err)

    def test_star_version_incompatible(self):
        rc, _, err = self.check(self.make_ref(version="2.5.0a"))
        self.assertEqual(rc, 1)
        self.assertIn("incompatible", err)

    def test_contig_mismatch(self):
        rc, _, err = self.check(self.make_ref(), contigs=("chr1", "chr2"))
        self.assertEqual(rc, 1)
        self.assertIn("contigs", err)

    def test_gene_not_in_gtf(self):
        rc, _, err = self.check(self.make_ref(), ann=("G1",))
        self.assertEqual(rc, 1)
        self.assertIn("not in the GTF", err)

    def test_extra_gtf_genes_only_warn(self):
        rc, _, err = self.check(self.make_ref(), ann=("G1", "G2", "G3"))
        self.assertEqual(rc, 0)
        self.assertIn("WARNING", err)

    def test_overhang_warning(self):
        rc, _, err = self.check(self.make_ref(overhang=49), rl=100)
        self.assertEqual(rc, 0)
        self.assertIn("sjdbOverhang 49", err)

    def test_no_overhang_warning_when_ok(self):
        rc, _, err = self.check(self.make_ref(overhang=99), rl=100)
        self.assertEqual((rc, err), (0, ""))

    def test_missing_dir_and_two_prefixes(self):
        self.assertEqual(self.check(self.p("nope"))[0], 1)
        d = self.make_ref()
        open(os.path.join(d, "other.grp"), "w").close()
        self.assertEqual(self.check(d)[0], 1)

    def test_usage(self):
        self.assertEqual(run("check_reference.py")[0], 2)


class TestHash(Tmp):
    def test_file_and_dir(self):
        f = self.write("a.txt", "hello\n")
        d = self.p("d")
        os.makedirs(d)
        open(os.path.join(d, "genomeParameters.txt"), "w").write("x")
        open(os.path.join(d, "huge.bin"), "w").write("y")
        rc, _, _ = run("hash_inputs.py", "--out", self.p("o.json"), "f=" + f, "r=" + d)
        self.assertEqual(rc, 0)
        j = json.loads(self.read("o.json"))
        self.assertEqual(j["f"]["sha256"],
                         "5891b5b522d5df086d0ff0b110fbd9d21bb4fc7163af34d08286a2e846f6be03")
        self.assertEqual(list(j["r"]["files"]), ["genomeParameters.txt"])

    def test_errors(self):
        self.assertEqual(run("hash_inputs.py")[0], 2)
        self.assertEqual(run("hash_inputs.py", "--out", self.p("o"), "x=" + self.p("nope"))[0], 1)


if __name__ == "__main__":
    unittest.main(verbosity=1)
