#!/usr/bin/env python3
"""Maximum read length over the first N reads of FASTQ file(s)."""
import gzip
import sys

VERSION = "0.1.0"
HELP = """\
Usage: max_read_length.py [--n_reads N] FASTQ [FASTQ ...]

Print the maximum read length found in the first N reads (default 10000) of every given
FASTQ file (.gz or plain). For paired-end data pass both mates: the maximum per-mate length
is reported. The pipeline uses (this value - 1) as STAR sjdbOverhang.

Options:
  --n_reads N   reads to inspect per file (default 10000)
  -h, --help    show this help and exit 0
  --version     print the version and exit 0

Exit codes: 0 ok, 1 unreadable/empty FASTQ, 2 usage error.

Example:
  max_read_length.py sample_R1.fq.gz sample_R2.fq.gz
"""


def longest(path, n):
    opener = gzip.open if path.endswith(".gz") else open
    best, reads = 0, 0
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                best = max(best, len(line.rstrip("\r\n")))
                reads += 1
                if reads >= n:
                    break
    return best, reads


def main(argv):
    n, files, i = 10000, [], 0
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(HELP)
            return 0
        if a == "--version":
            print("max_read_length.py " + VERSION)
            return 0
        if a == "--n_reads" and i + 1 < len(argv) and argv[i + 1].isdigit() and int(argv[i + 1]) > 0:
            n = int(argv[i + 1])
            i += 2
            continue
        if a.startswith("-"):
            sys.stderr.write("ERROR: bad argument: %s\n\n%s" % (a, HELP.split("\n\n")[0] + "\n"))
            return 2
        files.append(a)
        i += 1
    if not files:
        sys.stderr.write("ERROR: no FASTQ given\n\n%s" % (HELP.split("\n\n")[0] + "\n"))
        return 2
    best = 0
    for f in files:
        try:
            length, reads = longest(f, n)
        except (OSError, EOFError) as exc:
            sys.stderr.write("ERROR: cannot read %s: %s\n" % (f, exc))
            return 1
        if reads == 0:
            sys.stderr.write("ERROR: no reads in %s\n" % f)
            return 1
        best = max(best, length)
    print(best)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
