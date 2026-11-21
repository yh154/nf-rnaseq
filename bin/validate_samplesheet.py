#!/usr/bin/env python3
"""Validate the samplesheet, contrasts file and design formula (standard library only)."""
import csv
import os
import re
import sys
from fractions import Fraction

VERSION = "0.1.0"

HELP = """\
Usage: validate_samplesheet.py --samplesheet FILE --contrasts FILE [options]

Validate a bulk RNA-seq samplesheet, contrasts file and design formula before any heavy
process runs. All problems are collected and reported together.

Required arguments:
  --samplesheet FILE        CSV: sample,fastq_1,fastq_2,strandedness,condition[,batch,...]
  --contrasts FILE          either CSV id,variable,treatment,control (header) or a plain list of
                            pairs, one per line: treatment,control[,variable]  (comma or tab; '#'
                            comments; optional header, which is a first row made only of the words
                            id/variable/treatment/control; variable defaults to condition; the id is
                            <treatment>_vs_<control>). log2FC = treatment/control

Options:
  --design STR              DESeq2 design, additive terms only (default: "~ condition")
  --base_dir DIR            Directory that relative FASTQ paths are resolved against
                            (default: directory of the samplesheet)
  --no_file_check           Do not test that FASTQ files exist / match their extension
                            (used inside containers where the paths are not mounted)
  --out_samplesheet FILE    Write the normalised samplesheet (absolute paths)
  --out_contrasts FILE      Write the normalised contrasts (adds run_deseq2 column)
  --out_warnings FILE       Write warnings, one per line
  -h, --help                Show this help and exit 0
  --version                 Print the version and exit 0

Rules (each violation is an error):
  * required columns present; sample names match [A-Za-z0-9._-]+
  * a sample may occupy several rows (lanes); all its other columns must agree
  * single-end (empty fastq_2) and paired-end rows are not mixed
  * strandedness in {forward, reverse, none}
  * FASTQ: .fastq/.fq [.gz]; exist; unique paths and unique file names; same compression
    within a sample
  * design: additive terms of existing, non-empty columns; full rank model matrix with
    residual degrees of freedom (all design variables are treated as factors)
  * contrasts: unique ids; variable is in the design; treatment/control levels exist and
    differ; names match [A-Za-z0-9._-]+; no two contrasts give the same output name;
    at least one contrast remains after unreplicated samples are skipped (see below)

Unreplicated samples: a sample that is the only one with its level of a contrast variable
(e.g. the only sample with condition=drugC) is removed from the whole analysis (it is not
aligned, quantified, size-factor-normalised, tested, plotted or reported) with a warning,
repeated until every remaining level has >= 2 samples. A contrast that names a level that
was removed this way is dropped with a warning.

Exit codes: 0 valid, 1 validation failed, 2 usage error.

Example:
  validate_samplesheet.py --samplesheet samples.csv --contrasts contrasts.csv \\
      --design "~ batch + condition" --out_samplesheet samplesheet.valid.csv
"""

REQUIRED = ["sample", "fastq_1", "fastq_2", "strandedness", "condition"]
CONTRAST_COLS = ["id", "variable", "treatment", "control"]
STRAND = {"forward", "reverse", "none"}
NAME_RE = re.compile(r"^[A-Za-z0-9._-]+$")
FQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$")


class Usage(Exception):
    pass


def parse_args(argv):
    opts = {"design": "~ condition", "base_dir": None,
            "file_check": True, "samplesheet": None, "contrasts": None,
            "out_samplesheet": None, "out_contrasts": None, "out_warnings": None}
    flags = {"--no_file_check": ("file_check", False)}
    valued = {"--samplesheet", "--contrasts", "--design", "--base_dir",
              "--out_samplesheet", "--out_contrasts", "--out_warnings"}
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(HELP)
            sys.exit(0)
        if a == "--version":
            print("validate_samplesheet.py " + VERSION)
            sys.exit(0)
        if a in flags:
            key, val = flags[a]
            opts[key] = val
        elif a in valued:
            if i + 1 >= len(argv):
                raise Usage("missing value for " + a)
            opts[a[2:]] = argv[i + 1]
            i += 1
        else:
            raise Usage("unknown argument: " + a)
        i += 1
    for req in ("samplesheet", "contrasts"):
        if not opts[req]:
            raise Usage("--%s is required" % req)
    return opts


def read_csv(path, errors, what):
    try:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            rows = list(csv.reader(fh))
    except OSError as exc:
        errors.append("%s: cannot read %s (%s)" % (what, path, exc.strerror))
        return None, []
    rows = [[c.strip() for c in r] for r in rows if any(c.strip() for c in r)]
    if not rows:
        errors.append("%s is empty: %s" % (what, path))
        return None, []
    header = rows[0]
    if len(set(header)) != len(header):
        errors.append("%s: duplicated column names in header" % what)
    out = []
    for n, r in enumerate(rows[1:], start=2):
        if len(r) != len(header):
            errors.append("%s line %d: %d fields, header has %d" % (what, n, len(r), len(header)))
            continue
        d = dict(zip(header, r))
        d["_line"] = n
        out.append(d)
    return header, out


def is_gz(path):
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def rank(matrix):
    """Exact rank of a list of rows (Fractions)."""
    m = [[Fraction(x) for x in row] for row in matrix]
    r = 0
    ncol = len(m[0]) if m else 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c] / m[r][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        r += 1
    return r


def parse_design(design, errors):
    d = design.strip()
    if not d.startswith("~"):
        errors.append('design must start with "~" (got "%s")' % design)
        return []
    terms = [t.strip() for t in d[1:].split("+")]
    bad = [t for t in terms if not re.match(r"^[A-Za-z][A-Za-z0-9._]*$", t)]
    if not terms or bad:
        errors.append('design "%s": only additive terms of column names are supported '
                      "(e.g. \"~ batch + condition\")" % design)
        return []
    if len(set(terms)) != len(terms):
        errors.append("design repeats a term")
    return terms


def read_contrasts(path, errors):
    """Contrasts file: full format (id,variable,treatment,control with a header) or a plain list of
    pairs "treatment,control[,variable]" (comma or tab separated, optional header, '#' comments).
    Returns (CONTRAST_COLS, rows) with id = <treatment>_vs_<control> and variable = condition
    filled in where the file does not give them."""
    try:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            raw = fh.read().splitlines()
    except OSError as exc:
        errors.append("contrasts: cannot read %s (%s)" % (path, exc.strerror))
        return None, []
    lines = []
    for n, line in enumerate(raw, start=1):
        line = line.split("#", 1)[0].strip()
        if line:
            lines.append((n, line))
    if not lines:
        errors.append("contrasts file is empty: %s" % path)
        return None, []
    delim = "\t" if "\t" in lines[0][1] else ","
    table = [(n, [c.strip() for c in next(csv.reader([l], delimiter=delim))]) for n, l in lines]
    first = [c.lower() for c in table[0][1]]
    names = ("id", "variable", "treatment", "control")
    if all(c in names for c in first):                       # a header row: every cell is a column name
        if len(set(first)) != len(first):
            errors.append("contrasts: duplicated column names in header")
        missing = [c for c in ("treatment", "control") if c not in first]
        if missing:
            errors.append("contrasts: missing required column(s): " + ", ".join(missing))
            return CONTRAST_COLS, []
        cols, body = first, table[1:]
    else:                                                    # plain pairs
        cols, body = ["treatment", "control", "variable"][:max(2, len(table[0][1]))], table
    rows = []
    for n, f in body:
        if len(f) != len(cols):
            errors.append("contrasts line %d: %d fields, expected %d (%s)" % (n, len(f), len(cols), ",".join(cols)))
            continue
        d = dict(zip(cols, f))
        d.setdefault("variable", "")
        d["variable"] = d["variable"] or "condition"
        d["id"] = d.get("id") or "%s_vs_%s" % (d["treatment"], d["control"])
        d["_line"] = n
        rows.append(d)
    return CONTRAST_COLS, rows


def validate(opts):
    errors, warnings = [], []
    header, rows = read_csv(opts["samplesheet"], errors, "samplesheet")
    c_header, c_rows = read_contrasts(opts["contrasts"], errors)
    base = opts["base_dir"] or os.path.dirname(os.path.abspath(opts["samplesheet"]))
    samples = {}  # name -> list of rows

    if header is not None:
        missing = [c for c in REQUIRED if c not in header]
        if missing:
            errors.append("samplesheet: missing required column(s): " + ", ".join(missing))
            rows = []

    # ---- samplesheet rows
    paired = set()
    seen_files = {}
    seen_names = {}
    for r in rows:
        where = "samplesheet line %d" % r["_line"]
        s = r["sample"]
        if not s:
            errors.append(where + ": empty sample name")
            continue
        if not NAME_RE.match(s):
            errors.append(where + ': sample name "%s" must match [A-Za-z0-9._-]+ '
                          "(no whitespace or special characters)" % s)
        if r["strandedness"] not in STRAND:
            errors.append(where + ': strandedness "%s" not in {forward, reverse, none}'
                          % r["strandedness"])
        if not r["condition"]:
            errors.append(where + ": empty condition")
        if not r["fastq_1"]:
            errors.append(where + ": empty fastq_1")
            continue
        paired.add(bool(r["fastq_2"]))
        for col in ("fastq_1", "fastq_2"):
            p = r[col]
            if not p:
                continue
            if not FQ_RE.search(p):
                errors.append(where + ": %s must end in .fastq/.fq[.gz]: %s" % (col, p))
            full = p if os.path.isabs(p) else os.path.join(base, p)
            full = os.path.normpath(full)
            r[col] = full
            if full in seen_files:
                errors.append(where + ": %s already used on line %d: %s"
                              % (col, seen_files[full], p))
            seen_files[full] = r["_line"]
            base_name = os.path.basename(full)
            if base_name in seen_names and seen_names[base_name] != full:
                errors.append(where + ": FASTQ file name %s is also used by %s; file names "
                              "must be unique (they name the FastQC reports)"
                              % (base_name, seen_names[base_name]))
            seen_names.setdefault(base_name, full)
            if opts["file_check"]:
                if not os.path.isfile(full):
                    errors.append(where + ": %s not found: %s" % (col, full))
                elif FQ_RE.search(p) and is_gz(full) != bool(re.search(r"\.gz$", p)):
                    errors.append(where + ": %s extension disagrees with content "
                                  "(gzip magic): %s" % (col, p))
        if r["fastq_1"] and r["fastq_1"] == r["fastq_2"]:
            errors.append(where + ": fastq_1 and fastq_2 are the same file")
        samples.setdefault(s, []).append(r)

    if len(paired) > 1:
        errors.append("samplesheet mixes single-end and paired-end rows; "
                      "run them as separate analyses")

    other_cols = [c for c in (header or []) if c not in ("sample", "fastq_1", "fastq_2")]
    for s, srows in samples.items():
        for c in other_cols:
            vals = {r.get(c, "") for r in srows}
            if len(vals) > 1:
                errors.append('sample "%s": column %s differs between its rows (%s)'
                              % (s, c, ", ".join(sorted(map(repr, vals)))))
        comp = {bool(re.search(r"\.gz$", r[c])) for r in srows
                for c in ("fastq_1", "fastq_2") if r[c]}
        if len(comp) > 1:
            errors.append('sample "%s": mixes gzip and plain FASTQ files' % s)

    # ---- skip unreplicated samples (only level of a contrast variable): removed from everything
    cvars = sorted({c["variable"] for c in c_rows if header is not None and c["variable"] in header})
    gone_levels = set()
    while True:
        drop = {}
        for v in cvars:
            n = {}
            for s, srows in samples.items():
                n[srows[0].get(v, "")] = n.get(srows[0].get(v, ""), 0) + 1
            for s, srows in samples.items():
                lv = srows[0].get(v, "")
                if lv and n[lv] < 2:
                    drop.setdefault(s, (v, lv))
        if not drop:
            break
        for s, (v, lv) in drop.items():
            warnings.append('sample "%s" skipped from the whole analysis: it is the only sample with '
                            '%s=%s (no replicates)' % (s, v, lv))
            gone_levels.add((v, lv))
            del samples[s]
    rows = [r for r in rows if r["sample"] in samples]

    if not samples and not errors:
        errors.append("no samples left after skipping unreplicated samples")

    # ---- design
    terms = parse_design(opts["design"], errors)
    meta = {}  # sample -> first row
    for s, srows in samples.items():
        meta[s] = srows[0]
    for t in terms:
        if header is not None and t not in header:
            errors.append('design variable "%s" is not a samplesheet column' % t)
        elif t in ("sample", "fastq_1", "fastq_2", "strandedness"):
            errors.append('design variable "%s" cannot be used' % t)
        elif any(not meta[s].get(t, "") for s in meta):
            errors.append('design variable "%s" has empty values' % t)
    design_ok = bool(terms) and not any(
        t not in (header or []) or any(not meta[s].get(t, "") for s in meta) for t in terms)
    if design_ok and meta:
        names = list(meta)
        cols = [[1] * len(names)]
        for t in terms:
            levels = sorted({meta[s][t] for s in names})
            if len(levels) < 2:
                errors.append('design variable "%s" has a single level' % t)
            for lv in levels[1:]:
                cols.append([1 if meta[s][t] == lv else 0 for s in names])
        matrix = [list(row) for row in zip(*cols)]
        rk = rank(matrix)
        if rk < len(cols):
            errors.append("design is not full rank (%d columns, rank %d): variables are "
                          "confounded or a level combination is missing" % (len(cols), rk))
        elif len(names) <= len(cols):
            errors.append("design has no residual degrees of freedom "
                          "(%d samples, %d coefficients)" % (len(names), len(cols)))

    # ---- contrasts
    out_names = {}
    ids = set()
    valid_contrasts = []
    for c in c_rows:
        where = "contrasts line %d" % c["_line"]
        ok = True
        if not c["id"] or not NAME_RE.match(c["id"]):
            errors.append(where + ': id "%s" must match [A-Za-z0-9._-]+' % c["id"])
            ok = False
        elif c["id"] in ids:
            errors.append(where + ': duplicated id "%s"' % c["id"])
            ok = False
        ids.add(c["id"])
        for col in ("treatment", "control"):
            if not c[col] or not NAME_RE.match(c[col]):
                errors.append(where + ': %s "%s" must match [A-Za-z0-9._-]+' % (col, c[col]))
                ok = False
        if c["treatment"] and c["treatment"] == c["control"]:
            errors.append(where + ": treatment and control are identical")
            ok = False
        if c["variable"] not in terms:
            errors.append(where + ': variable "%s" is not in the design (%s)'
                          % (c["variable"], opts["design"]))
            ok = False
        elif meta and header is not None and c["variable"] in header:
            counts = {}
            for s in meta:
                counts.setdefault(meta[s][c["variable"]], []).append(s)
            for col in ("treatment", "control"):
                lv = c[col]
                if lv not in counts and (c["variable"], lv) in gone_levels:
                    warnings.append(where + ': contrast "%s" skipped: level "%s" of %s has no '
                                    "replicated samples" % (c["id"], lv, c["variable"]))
                    ok = False
                elif lv not in counts:
                    errors.append(where + ': %s level "%s" not found in column %s '
                                  "(levels: %s)" % (col, lv, c["variable"],
                                                    ", ".join(sorted(counts))))
                    ok = False
            c["run_deseq2"] = "true"
        name = "%s_vs_%s" % (c["treatment"], c["control"])
        if name in out_names:
            errors.append(where + ': output name "%s" duplicates line %d'
                          % (name, out_names[name]))
            ok = False
        out_names[name] = c["_line"]
        if ok:
            valid_contrasts.append(c)
    if not c_rows and c_header is not None and not errors:
        errors.append("contrasts file has no rows")
    elif c_rows and not valid_contrasts and not errors:
        errors.append("no contrast is left after skipping unreplicated samples (see warnings)")

    return errors, warnings, header, rows, c_header, valid_contrasts


def write_outputs(opts, errors, warnings, header, rows, c_header, contrasts):
    if opts["out_warnings"]:
        with open(opts["out_warnings"], "w") as fh:
            fh.write("".join(w + "\n" for w in warnings))
    if errors:
        return
    if opts["out_samplesheet"]:
        with open(opts["out_samplesheet"], "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            for r in rows:
                w.writerow([r[c] for c in header])
    if opts["out_contrasts"]:
        cols = list(c_header) + ["run_deseq2"]
        with open(opts["out_contrasts"], "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for c in contrasts:
                w.writerow([c.get(k, "true" if k == "run_deseq2" else "") for k in cols])


def main(argv):
    try:
        opts = parse_args(argv)
    except Usage as exc:
        sys.stderr.write("ERROR: %s\n\n%s" % (exc, HELP.split("\n\n")[0] + "\n"))
        return 2
    errors, warnings, header, rows, c_header, contrasts = validate(opts)
    write_outputs(opts, errors, warnings, header, rows, c_header, contrasts)
    for w in warnings:
        sys.stderr.write("WARNING: %s\n" % w)
    if errors:
        for e in errors:
            sys.stderr.write("ERROR: %s\n" % e)
        sys.stderr.write("Validation failed: %d error(s)\n" % len(errors))
        return 1
    nsamp = len({r["sample"] for r in rows})
    layout = "paired-end" if rows and rows[0]["fastq_2"] else "single-end"
    print("Validation OK: %d samples (%d FASTQ rows), %s, %d contrast(s)"
          % (nsamp, len(rows), layout, len(contrasts)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
