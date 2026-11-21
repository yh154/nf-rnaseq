#!/usr/bin/env python3
"""Unit tests for bin/validate_samplesheet.py (valid case and every failure)."""
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "bin", "validate_samplesheet.py")
DATA = os.path.join(ROOT, "tests", "data")


def run(*args):
    p = subprocess.run([sys.executable, "-I", SCRIPT, *args], capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.sheet = os.path.join(DATA, "samples_pe.csv")
        self.contr = os.path.join(DATA, "contrasts.csv")

    def edit_sheet(self, fn, name="s.csv"):
        with open(self.sheet, newline="") as fh:
            rows = list(csv.DictReader(fh))
            cols = list(rows[0])
        rows = fn(rows) or rows
        cols = list(rows[0]) if rows else cols
        path = os.path.join(self.tmp, name)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        return path

    def write_contrasts(self, text):
        path = os.path.join(self.tmp, "c.csv")
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def check(self, sheet=None, contr=None, extra=(), ok=False, msg=None):
        rc, out, err = run("--samplesheet", sheet or self.sheet, "--contrasts", contr or self.contr,
                           "--base_dir", DATA, *extra)
        if ok:
            self.assertEqual(rc, 0, err)
        else:
            self.assertEqual(rc, 1, out + err)
            if msg:
                self.assertIn(msg, err)
        return rc, out, err


class TestValid(Base):
    def test_pe_valid_and_outputs(self):
        s, c, w = (os.path.join(self.tmp, n) for n in ("s.csv", "c.csv", "w.txt"))
        self.check(ok=True, extra=("--design", "~ batch + condition", "--out_samplesheet", s,
                                   "--out_contrasts", c, "--out_warnings", w))
        rows = list(csv.DictReader(open(s)))
        self.assertTrue(all(os.path.isabs(r["fastq_1"]) for r in rows))
        self.assertEqual(open(c).read().splitlines()[1],
                         "treated_vs_control,condition,treated,control,true")
        self.assertEqual(open(w).read(), "")

    def test_se_valid(self):
        self.check(sheet=os.path.join(DATA, "samples_se.csv"), ok=True)

    def test_lanes_allowed(self):
        _, out, _ = self.check(ok=True)
        self.assertIn("6 samples (7 FASTQ rows)", out)


class TestSamplesheetFailures(Base):
    def test_missing_file(self):
        def f(rows): rows[0]["fastq_2"] = "fastq/nope.fq.gz"
        self.check(self.edit_sheet(f), msg="not found")

    def test_no_file_check_skips_existence(self):
        def f(rows): rows[0]["fastq_2"] = "fastq/nope.fq.gz"
        self.check(self.edit_sheet(f), extra=("--no_file_check",), ok=True)

    def test_mixed_se_pe(self):
        def f(rows): rows[0]["fastq_2"] = ""
        self.check(self.edit_sheet(f), msg="mixes single-end and paired-end")

    def test_bad_strandedness(self):
        def f(rows): rows[0]["strandedness"] = "unstranded"
        self.check(self.edit_sheet(f), msg="strandedness")

    def test_bad_sample_name(self):
        def f(rows):
            for r in rows:
                if r["sample"] == "control_rep1":
                    r["sample"] = "control rep1"
        self.check(self.edit_sheet(f), msg="must match")

    def test_lane_rows_disagree(self):
        def f(rows):
            for r in rows:
                if r["fastq_1"].endswith("L002_pe_R1.fq.gz"):
                    r["condition"] = "treated"
        self.check(self.edit_sheet(f), msg="differs between its rows")

    def test_duplicate_fastq(self):
        def f(rows): rows[1]["fastq_1"] = rows[0]["fastq_1"]
        self.check(self.edit_sheet(f), msg="already used")

    def test_duplicate_basename_in_different_dirs(self):
        os.makedirs(os.path.join(self.tmp, "other"))
        shutil.copy(os.path.join(DATA, "fastq", "treated_rep1_pe_R1.fq.gz"),
                    os.path.join(self.tmp, "other", "control_rep1_pe_R1.fq.gz"))
        def f(rows): rows[4]["fastq_1"] = os.path.join(self.tmp, "other", "control_rep1_pe_R1.fq.gz")
        self.check(self.edit_sheet(f), msg="must be unique")

    def test_wrong_extension(self):
        def f(rows): rows[0]["fastq_1"] = "genome.fa"
        self.check(self.edit_sheet(f), msg=".fastq/.fq")

    def test_gzip_magic_mismatch(self):
        shutil.copy(os.path.join(DATA, "genome.fa"), os.path.join(self.tmp, "fake.fq.gz"))
        def f(rows): rows[0]["fastq_1"] = os.path.join(self.tmp, "fake.fq.gz")
        self.check(self.edit_sheet(f), msg="gzip magic")

    def test_missing_required_column(self):
        def f(rows):
            for r in rows:
                del r["strandedness"]
        self.check(self.edit_sheet(f), msg="missing required column")

    def test_empty_file(self):
        p = os.path.join(self.tmp, "empty.csv")
        open(p, "w").close()
        self.check(p, msg="empty")


class TestDesignAndContrasts(Base):
    def test_all_but_one_level_unreplicated_fails(self):
        def f(rows): return [r for r in rows if r["sample"] not in ("treated_rep2", "treated_rep3")]
        _, _, err = self.check(self.edit_sheet(f), extra=("--design", "~ condition"),
                               msg="has a single level")
        self.assertIn('sample "treated_rep1" skipped', err)

    def test_unreplicated_sample_is_skipped_everywhere(self):
        def f(rows):
            return rows + [dict(rows[0], sample="lone", condition="other", fastq_1="fastq/lone_R1.fq.gz",
                                fastq_2="fastq/lone_R2.fq.gz")]
        s, w = os.path.join(self.tmp, "s.csv"), os.path.join(self.tmp, "w.txt")
        _, out, err = self.check(self.edit_sheet(f), ok=True, extra=(
            "--no_file_check", "--out_samplesheet", s, "--out_warnings", w))
        self.assertIn('sample "lone" skipped from the whole analysis', err)
        self.assertIn("6 samples (7 FASTQ rows)", out)
        self.assertNotIn("lone", {r["sample"] for r in csv.DictReader(open(s))})
        self.assertIn("lone", open(w).read())

    def test_contrast_with_skipped_level_is_dropped(self):
        def f(rows):
            return rows + [dict(rows[0], sample="lone", condition="other", fastq_1="fastq/lone_R1.fq.gz",
                                fastq_2="fastq/lone_R2.fq.gz")]
        c = self.write_contrasts("treated,control\nother,control\n")
        o = os.path.join(self.tmp, "oc.csv")
        _, _, err = self.check(self.edit_sheet(f), c, ok=True,
                               extra=("--no_file_check", "--out_contrasts", o))
        self.assertIn('contrast "other_vs_control" skipped', err)
        self.assertEqual(len(open(o).read().splitlines()), 2)

    def test_confounded_design(self):
        def f(rows):
            for r in rows:
                r["batch"] = "x" if r["condition"] == "treated" else "y"
        self.check(self.edit_sheet(f), extra=("--design", "~ batch + condition"),
                   msg="not full rank")

    def test_design_variable_missing(self):
        self.check(extra=("--design", "~ nonexistent + condition"), msg="not a samplesheet column")

    def test_design_unsupported_syntax(self):
        self.check(extra=("--design", "~ batch * condition"), msg="only additive terms")

    def test_no_degrees_of_freedom(self):
        def f(rows):
            for r in rows:
                r["batch"] = r["sample"]
        self.check(self.edit_sheet(f), extra=("--design", "~ batch + condition"))

    def test_contrast_variable_not_in_design(self):
        c = self.write_contrasts("id,variable,treatment,control\nb,batch,b2,b1\n")
        self.check(contr=c, extra=("--design", "~ condition"), msg="not in the design")

    def test_contrast_batch_in_design_ok(self):
        c = self.write_contrasts("id,variable,treatment,control\nb,batch,b2,b1\n"
                                 "t,condition,treated,control\n")
        self.check(contr=c, extra=("--design", "~ batch + condition"), ok=True)

    def test_contrast_level_missing(self):
        c = self.write_contrasts("id,variable,treatment,control\nx,condition,drug,control\n")
        self.check(contr=c, msg='level "drug" not found')

    def test_contrast_same_levels(self):
        c = self.write_contrasts("id,variable,treatment,control\nx,condition,control,control\n")
        self.check(contr=c, msg="identical")

    def test_swapped_contrast_changes_output_name_not_validity(self):
        c = self.write_contrasts("id,variable,treatment,control\nx,condition,control,treated\n")
        out = os.path.join(self.tmp, "o.csv")
        self.check(contr=c, ok=True, extra=("--out_contrasts", out))
        self.assertIn("x,condition,control,treated", open(out).read())

    def test_duplicate_output_name(self):
        c = self.write_contrasts("id,variable,treatment,control\na,condition,treated,control\n"
                                 "b,condition,treated,control\n")
        self.check(contr=c, msg="duplicates")

    def test_duplicate_id(self):
        c = self.write_contrasts("id,variable,treatment,control\na,condition,treated,control\n"
                                 "a,condition,control,treated\n")
        self.check(contr=c, msg="duplicated id")

    def test_contrasts_missing_column(self):
        c = self.write_contrasts("id,variable,treatment\na,condition,treated\n")
        self.check(contr=c, msg="missing required column")


class TestContrastFormats(Base):
    """--contrast accepts a plain list of pairs as well as the full CSV."""

    def pairs(self, text):
        return self.write_contrasts(text)

    def out(self, contr, extra=()):
        o = os.path.join(self.tmp, "o.csv")
        self.check(contr=contr, ok=True, extra=("--out_contrasts", o, *extra))
        with open(o) as fh:
            return fh.read().splitlines()[1:]

    def test_plain_pairs_comma_tab_comments(self):
        for text in ("treated,control\n", "treated\tcontrol\n", "# note\n\ntreated , control  # trailing\n"):
            self.assertEqual(self.out(self.pairs(text)), ["treated_vs_control,condition,treated,control,true"])

    def test_pairs_with_header_and_variable_column(self):
        self.assertEqual(self.out(self.pairs("treatment,control\ntreated,control\n")),
                         ["treated_vs_control,condition,treated,control,true"])
        lines = self.out(self.pairs("b2,b1,batch\ntreated,control,condition\n"),
                         extra=("--design", "~ batch + condition"))
        self.assertEqual(lines, ["b2_vs_b1,batch,b2,b1,true", "treated_vs_control,condition,treated,control,true"])

    def test_level_named_control_is_not_a_header(self):
        # "treated,control" contains the word control but is data, not a header
        self.assertEqual(len(self.out(self.pairs("treated,control\n"))), 1)

    def test_swapped_pair_gives_opposite_direction_and_name(self):
        self.assertEqual(self.out(self.pairs("control,treated\n")), ["control_vs_treated,condition,control,treated,true"])

    def test_several_pairs_keep_order(self):
        c = self.pairs("treated,control\ncontrol,treated\n")
        self.assertEqual([l.split(",")[0] for l in self.out(c)], ["treated_vs_control", "control_vs_treated"])

    def test_full_format_still_works_with_comments(self):
        self.assertEqual(self.out(self.pairs("# c\nid,variable,treatment,control\nmine,condition,treated,control\n")),
                         ["mine,condition,treated,control,true"])

    def test_pair_errors(self):
        self.check(contr=self.pairs("treated\n"), msg="1 fields, expected 2")
        self.check(contr=self.pairs("treated,control\ntreated,control,condition,x\n"), msg="fields, expected")
        self.check(contr=self.pairs("treated,drug\n"), msg='control level "drug" not found')
        self.check(contr=self.pairs("treated,control\ntreated,control\n"), msg="duplicate")
        self.check(contr=self.pairs("treated,treated\n"), msg="identical")
        self.check(contr=self.pairs("# only a comment\n"), msg="empty")
        self.check(contr=self.pairs("id,variable,treatment\na,condition,treated\n"), msg="missing required column")
        self.check(contr=self.pairs("treated,control,nonexistent\n"), msg="not in the design")


class TestCli(Base):
    def test_help_and_version(self):
        rc, out, _ = run("--help")
        self.assertEqual(rc, 0)
        self.assertIn("Usage", out)
        rc, out, _ = run("--version")
        self.assertEqual(rc, 0)

    def test_usage_errors_exit_2(self):
        self.assertEqual(run()[0], 2)
        self.assertEqual(run("--bogus")[0], 2)
        self.assertEqual(run("--samplesheet")[0], 2)


if __name__ == "__main__":
    unittest.main(verbosity=1)
