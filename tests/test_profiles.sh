#!/usr/bin/env bash
# The execution profiles resolve to the intended configuration (no cluster or container runtime
# needed): slurm + apptainer, local + docker, and the pinned image references.
# Usage: tests/test_profiles.sh [nextflow]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
NF="${1:-nextflow}"
fail() { echo "FAIL: $*" >&2; exit 1; }

hpc="$($NF config -profile slurm,apptainer . 2>/dev/null)"
grep -q "executor = 'slurm'" <<< "$hpc"       || fail "slurm profile: executor is not slurm"
# --queue cannot be given to `nextflow config`: evaluate conf/slurm.config with params.queue set
tmp="$(mktemp -d "${TMPDIR:-/tmp}/profiles.XXXXXX")"
printf "params.queue = 'testq'\nincludeConfig '%s/conf/slurm.config'\n" "$ROOT" > "$tmp/q.config"
$NF -c "$tmp/q.config" config . 2>/dev/null | grep -q "queue = 'testq'" || fail "--queue is not applied to process.queue"
rm -rf "$tmp"
grep -q "enabled = true" <<< "$(sed -n '/^apptainer {/,/^}/p' <<< "$hpc")" || fail "apptainer not enabled"
grep -q "autoMounts = true" <<< "$hpc"         || fail "apptainer autoMounts not set"
grep -q "resourceLimits" <<< "$hpc"            || fail "no process.resourceLimits in the slurm profile"
grep -q "errorStrategy" <<< "$hpc"             || fail "no retry strategy"

loc="$($NF config -profile local,docker . 2>/dev/null)"
grep -q "executor = 'local'" <<< "$loc"        || fail "local profile: executor is not local"
grep -q "enabled = true" <<< "$(sed -n '/^docker {/,/^}/p' <<< "$loc")" || fail "docker not enabled"

# the same six references for both container profiles, digest-pinned, no placeholders
for cfg in "$hpc" "$loc"; do
    for label in fastqc star_rsem python deseq2 fgsea gsva; do
        grep -Eq "withLabel:$label .*container = 'quay.io/biocontainers/[a-z0-9._-]+@sha256:[0-9a-f]{64}'|container = 'quay.io/biocontainers/[a-z0-9._-]+@sha256:[0-9a-f]{64}'" <<< "$cfg" || fail "no digest for $label"
    done
    [[ "$(grep -c "container = 'quay.io/biocontainers/.*@sha256:[0-9a-f]\{64\}'" <<< "$cfg")" -eq 6 ]] || fail "expected 6 digest-pinned images"
done
[[ "$(grep -o "container = '[^']*'" <<< "$hpc" | sort)" == "$(grep -o "container = '[^']*'" <<< "$loc" | sort)" ]] \
    || fail "docker and apptainer profiles use different images"
if grep -rqiE "REPLACE_WITH|placeholder" conf/containers.config; then fail "placeholder in conf/containers.config"; fi
echo "profiles resolve: slurm+apptainer and local+docker, 6 identical digest-pinned images"
