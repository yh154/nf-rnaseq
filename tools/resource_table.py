#!/usr/bin/env python3
"""Markdown resource table from a Nextflow trace.txt (standard library only)."""
import csv
import re
import sys

VERSION = "0.1.0"
HELP = """\
Usage: resource_table.py TRACE.txt

Summarise a Nextflow trace per process: tasks, requested CPUs, peak RSS and wall time (the
maximum over the tasks of that process). Prints a Markdown table; used for the README.

Options: -h, --help / --version

Exit codes: 0 ok, 1 unreadable trace, 2 usage error.
Example: resource_table.py results/pipeline_info/trace.txt
"""


def to_mb(s):
    m = re.match(r"([\d.]+)\s*(KB|MB|GB|TB|B)", s or "")
    if not m:
        return 0.0
    return float(m.group(1)) * {"B": 1 / 1e6, "KB": 1 / 1e3, "MB": 1, "GB": 1e3, "TB": 1e6}[m.group(2)]


def to_s(s):
    t = 0.0
    for v, u in re.findall(r"([\d.]+)\s*(ms|s|m|h|d)", s or ""):
        t += float(v) * {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}[u]
    return t


def main(argv):
    if len(argv) == 1 and argv[0] in ("-h", "--help"):
        print(HELP)
        return 0
    if len(argv) == 1 and argv[0] == "--version":
        print("resource_table.py " + VERSION)
        return 0
    if len(argv) != 1:
        sys.stderr.write("Usage: resource_table.py TRACE.txt\n")
        return 2
    try:
        with open(argv[0], newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    except OSError as exc:
        sys.stderr.write("ERROR: %s\n" % exc)
        return 1
    agg = {}
    for r in rows:
        name = r["name"].split(" (")[0]
        a = agg.setdefault(name, {"n": 0, "cpus": r.get("cpus", ""), "mb": 0.0, "s": 0.0})
        a["n"] += 1
        a["mb"] = max(a["mb"], to_mb(r.get("peak_rss")))
        a["s"] = max(a["s"], to_s(r.get("realtime")))
    print("| process | tasks | CPUs | peak RSS | wall time (max) |")
    print("|---|---|---|---|---|")
    for name, a in agg.items():
        rss = "%.1f GB" % (a["mb"] / 1000) if a["mb"] >= 1000 else "%.0f MB" % a["mb"]
        print("| %s | %d | %s | %s | %ds |" % (name, a["n"], a["cpus"], rss, round(a["s"])))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
