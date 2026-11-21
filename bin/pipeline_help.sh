#!/usr/bin/env bash
# Print the exact local and HPC commands and resource requirements.
set -euo pipefail
VERSION="0.1.0"
usage() {
cat <<'USAGE'
Usage: pipeline_help.sh [-h|--help|--version]

Print the exact commands for a Slurm + Apptainer run (primary target) and a local Docker
run (development/tests), and the default per-process resource requirements.

Slurm + Apptainer:
  nextflow run <pipeline_dir> -profile slurm,apptainer --queue <partition> -c site.config \
      --input samples.csv --contrast contrasts.csv --gfasta genome.fa --gtf genome.gtf \
      --outdir results
Local Docker (tests, small data):
  nextflow run <pipeline_dir> -profile local,docker --input samples.csv ...

Resources: STAR_RSEM 8 CPU / 48 GB / 12 h (human); other steps are small.
Site caps: process.resourceLimits in your site config (see conf/site.example.config).
USAGE
}
case "${1:-}" in
  "") usage ;;
  -h|--help) usage ;;
  --version) echo "pipeline_help.sh $VERSION" ;;
  *) usage >&2; exit 2 ;;
esac
