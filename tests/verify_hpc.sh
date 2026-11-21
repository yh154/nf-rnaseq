#!/usr/bin/env bash
# Manual verification on a Slurm cluster with Apptainer (cannot run in CI or on a laptop).
# Runs the synthetic test dataset with -profile slurm,apptainer,test, then every output check.
#
# Usage: tests/verify_hpc.sh --queue PARTITION [--site-config FILE] [--outdir DIR]
#   needs on the login node: nextflow >= 25.04 (Java 17+), apptainer, python3, outbound access to
#   quay.io for the first image pull (or a pre-filled NXF_APPTAINER_CACHEDIR).
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 1
QUEUE=""; SITE=""; OUT="$ROOT/verify_hpc_out"
usage() { echo "Usage: tests/verify_hpc.sh --queue PARTITION [--site-config FILE] [--outdir DIR]"; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --queue) QUEUE="$2"; shift ;;
        --site-config) SITE="$2"; shift ;;
        --outdir) OUT="$2"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
    shift
done
[[ -n "$QUEUE" ]] || { usage >&2; exit 2; }
command -v apptainer >/dev/null || { echo "ERROR: apptainer not found" >&2; exit 1; }
command -v nextflow >/dev/null || { echo "ERROR: nextflow not found" >&2; exit 1; }

python3 -I tests/make_synthetic.py >/dev/null
echo "== image check (same digests as docker) =="
failed=0
for label in fastqc star_rsem python deseq2 fgsea gsva; do
    img="$(grep "withLabel: $label " conf/containers.config | sed "s/.*container = '\(.*\)'.*/\1/")"
    if apptainer exec "docker://$img" true >/dev/null 2>&1; then echo "ok   $label"; else echo "FAIL $label ($img)"; failed=1; fi
done

echo "== pipeline on the synthetic data =="
cfg=(); [[ -n "$SITE" ]] && cfg=(-c "$SITE")
nextflow run . -profile slurm,apptainer,test --queue "$QUEUE" "${cfg[@]}" --outdir "$OUT" -ansi-log false || failed=1

echo "== output checks =="
chk() { if "$@"; then :; else echo "FAIL: $*"; failed=1; fi; }
chk python3 -I tests/check_phase3.py "$OUT" tests/data/samples_pe.csv
chk python3 -I tests/check_phase4.py tables "$OUT/expression" "$OUT" tests/data/samples_pe.csv
chk python3 -I tests/check_phase4.py deseq "$OUT/diffexp" "$OUT/expression" tests/data/samples_pe.csv
chk python3 -I tests/check_phase4.py golden "$OUT/expression"
chk python3 -I tests/check_phase5.py "$OUT" tests/data/samples_pe.csv 1 "$(date +%d%m%y)" 2
chk python3 -I tests/check_phase6.py "$OUT/gsea" "$OUT/gsva"

echo
echo "### Slurm + Apptainer verification"
echo "- date: $(date +%Y-%m-%d), host: $(hostname)"
echo "- Nextflow: $(nextflow -version 2>/dev/null | grep -m1 version | xargs)"
echo "- Apptainer: $(apptainer --version)"
echo "- command: nextflow run . -profile slurm,apptainer,test --queue $QUEUE ${SITE:+-c $SITE }--outdir $OUT"
echo "- executor in trace: $(awk -F'\t' 'NR==2 {print $3}' "$OUT/pipeline_info/trace.txt" 2>/dev/null) (native_id = Slurm job id)"
[[ $failed -eq 0 ]] && echo "- outcome: PASSED" || echo "- outcome: FAILED (see output above)"
exit $failed
