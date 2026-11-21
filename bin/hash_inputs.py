#!/usr/bin/env python3
"""sha256 of the reference inputs, for provenance."""
import hashlib
import json
import os
import sys

VERSION = "0.1.0"
HELP = """\
Usage: hash_inputs.py --out FILE LABEL=PATH [LABEL=PATH ...]

Write a JSON file with size and sha256 for each input (provenance only, never compared).
A PATH that is a directory (a RSEM reference) is summarised by its small key files
(genomeParameters.txt, *.ti, *.transcripts.fa) rather than hashing the multi-GB STAR index.

Options:
  --out FILE   output JSON
  -h, --help   show this help and exit 0
  --version    print the version and exit 0

Example:
  hash_inputs.py --out inputs_sha256.json gfasta=genome.fa gtf=genome.gtf
"""


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def entry(path):
    p = os.path.realpath(path)
    if os.path.isdir(p):
        keys = sorted(f for f in os.listdir(p)
                      if f == "genomeParameters.txt" or f.endswith((".ti", ".transcripts.fa")))
        return {"directory": True,
                "files": {k: {"size": os.path.getsize(os.path.join(p, k)),
                              "sha256": sha(os.path.join(p, k))} for k in keys}}
    return {"size": os.path.getsize(p), "sha256": sha(p)}


def main(argv):
    out, items, i = None, [], 0
    while i < len(argv):
        a = argv[i]
        if a in ("-h", "--help"):
            print(HELP)
            return 0
        if a == "--version":
            print("hash_inputs.py " + VERSION)
            return 0
        if a == "--out" and i + 1 < len(argv):
            out = argv[i + 1]
            i += 2
            continue
        if "=" in a and not a.startswith("-"):
            items.append(a.split("=", 1))
            i += 1
            continue
        sys.stderr.write("ERROR: bad argument: %s\n\n%s" % (a, HELP.split("\n\n")[0] + "\n"))
        return 2
    if not out or not items:
        sys.stderr.write("ERROR: --out and at least one LABEL=PATH required\n\n%s"
                         % (HELP.split("\n\n")[0] + "\n"))
        return 2
    result = {}
    for label, path in items:
        if not os.path.exists(path):
            sys.stderr.write("ERROR: not found: %s\n" % path)
            return 1
        result[label] = entry(path)
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
