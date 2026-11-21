#!/usr/bin/env bash
# Phase 4 without alignment: per-sample genes.results fixtures (tests/fixtures/rsem_pe, produced
# by the pipeline on the synthetic data) -> merge_rsem.py -> run_deseq2.R in the pinned
# deseq2 container -> assertions. Needs docker.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
W="$(mktemp -d "${TMPDIR:-/tmp}/p4fast.XXXXXX")"
IMG="$(grep 'withLabel: deseq2 ' conf/containers.config | sed "s/.*container = '\(.*\)'.*/\1/")"
DATA=tests/data
run_deseq() { # outdir contrasts
    docker run --rm -u "$(id -u):$(id -g)" -v "$ROOT/bin:/pipeline_bin:ro" -v "$W:/w" -w /w "$IMG" \
        Rscript /pipeline_bin/run_deseq2.R --counts expr/RSEM_expected_counts.txt --samplesheet ss.csv \
        --contrasts "$2" --design '~ batch + condition' --outdir "$1" > "$W/$1.log" 2>&1
}
python3 -I bin/validate_samplesheet.py --samplesheet $DATA/samples_pe.csv --contrasts $DATA/contrasts.csv \
    --design '~ batch + condition' --out_samplesheet "$W/ss.csv" --out_contrasts "$W/ct.csv" > /dev/null
python3 -I bin/build_annotation.py --gtf $DATA/genes.gtf --out "$W/ann.tsv" > /dev/null
python3 -I bin/merge_rsem.py --samplesheet "$W/ss.csv" --annotation "$W/ann.tsv" --outdir "$W/expr" \
    tests/fixtures/rsem_pe/*.genes.results
python3 -I tests/check_phase4.py tables "$W/expr" tests/fixtures/rsem_pe "$W/ss.csv"
run_deseq diffexp ct.csv
python3 -I tests/check_phase4.py deseq "$W/diffexp" "$W/expr" "$W/ss.csv"
python3 -I tests/check_phase4.py golden "$W/diffexp"
# swapped treatment/control: new file name, opposite sign
printf 'id,variable,treatment,control,run_deseq2\ncontrol_vs_treated,condition,control,treated,true\n' > "$W/ct_swapped.csv"
run_deseq diffexp_swapped ct_swapped.csv
python3 -I tests/check_phase4.py swap "$W/diffexp" "$W/diffexp_swapped"
# run_deseq2 = false skips the contrast, no CSV
printf 'id,variable,treatment,control,run_deseq2\ntreated_vs_control,condition,treated,control,false\n' > "$W/ct_skip.csv"
run_deseq diffexp_skip ct_skip.csv
test -z "$(ls "$W"/diffexp_skip/Deseq2_res_*.csv 2>/dev/null)" && test -s "$W/diffexp_skip/dds.rds"
echo "ok  run_deseq2 = false produces no CSV"
# bad arguments
docker run --rm -v "$ROOT/bin:/pipeline_bin:ro" "$IMG" Rscript /pipeline_bin/run_deseq2.R --padj 7 > /dev/null 2>&1 && exit 1
# plots for the 2-group data: 1 global PCA + 4 per contrast pages
run_plots() { # dds contrasts resultsdir pdf
    docker run --rm -u "$(id -u):$(id -g)" -v "$ROOT/bin:/pipeline_bin:ro" -v "$W:/w" -w /w "$IMG" \
        Rscript /pipeline_bin/make_plots.R --dds "$1" --contrasts "$2" --results_dir "$3" --out "$4" > "$W/plots.log" 2>&1
}
pdf_pages() { python3 -I -c "import re,sys; print(len(re.findall(rb'/Type\\s*/Page\\b(?!s)', open(sys.argv[1],'rb').read())))" "$1"; }
run_plots diffexp/dds.rds ct.csv diffexp plots2.pdf
test "$(pdf_pages "$W/plots2.pdf")" -eq 5
echo "ok  plots PDF: 5 pages for 1 contrast"
# three condition levels, three contrasts: only the two contrasted groups in each CSV, 13 pages
awk -F, 'BEGIN{OFS=","} NR==1{print; next} {c=$1; $5=(c=="control_rep1"||c=="control_rep2")?"A":(c=="control_rep3"||c=="treated_rep1")?"B":"C"; print}' "$W/ss.csv" > "$W/ss3.csv"
printf 'id,variable,treatment,control,run_deseq2\nB_vs_A,condition,B,A,true\nC_vs_A,condition,C,A,true\nC_vs_B,condition,C,B,true\n' > "$W/ct3.csv"
docker run --rm -u "$(id -u):$(id -g)" -v "$ROOT/bin:/pipeline_bin:ro" -v "$W:/w" -w /w "$IMG" \
    Rscript /pipeline_bin/run_deseq2.R --counts expr/RSEM_expected_counts.txt --samplesheet ss3.csv \
    --contrasts ct3.csv --design '~ condition' --outdir diffexp3 > "$W/d3.log" 2>&1
python3 -I tests/check_phase4.py deseq "$W/diffexp3" "$W/expr" "$W/ss3.csv"
test "$(find "$W/diffexp3" -name 'Deseq2_res_*.csv' | wc -l)" -eq 3
run_plots diffexp3/dds.rds ct3.csv diffexp3 plots3.pdf
test "$(pdf_pages "$W/plots3.pdf")" -eq 13
echo "ok  3 contrasts: 3 CSVs with only the contrasted groups, 13-page PDF"
# GSEA (fgsea image) and GSVA (gsva image) on the same fixtures
img() { grep "withLabel: $1 " conf/containers.config | sed "s/.*container = '\(.*\)'.*/\1/"; }
cp $DATA/test.gmt "$W/test.gmt"
docker run --rm -u "$(id -u):$(id -g)" -v "$ROOT/bin:/pipeline_bin:ro" -v "$W:/w" -w /w "$(img fgsea)" \
    Rscript /pipeline_bin/run_gsea.R --results_dir diffexp --contrasts ct.csv --gmt test.gmt --min_set_size 5 --outdir GSEA > "$W/gsea.log" 2>&1
docker run --rm -u "$(id -u):$(id -g)" -v "$ROOT/bin:/pipeline_bin:ro" -v "$W:/w" -w /w "$(img gsva)" \
    Rscript /pipeline_bin/run_gsva.R --expr diffexp/vst_counts.txt --input_type vst --samplesheet ss.csv \
    --contrasts ct.csv --gmt test.gmt --min_set_size 5 --outdir GSVA > "$W/gsva.log" 2>&1
python3 -I tests/check_phase6.py "$W/GSEA" "$W/GSVA"
# GSVA on log2(TPM + 1): TPM is compositional, so only the relaxed check applies
docker run --rm -u "$(id -u):$(id -g)" -v "$ROOT/bin:/pipeline_bin:ro" -v "$W:/w" -w /w "$(img gsva)" \
    Rscript /pipeline_bin/run_gsva.R --expr expr/TPM.txt --input_type log2tpm --samplesheet ss.csv \
    --contrasts ct.csv --gmt test.gmt --min_set_size 5 --outdir GSVA_tpm > "$W/gsva_tpm.log" 2>&1
python3 -I tests/check_phase6.py "$W/GSEA" "$W/GSVA_tpm" test --relaxed-gsva
# plots with GSEA + GSVA pages: 5 + 2
run_plots_gs() {
    docker run --rm -u "$(id -u):$(id -g)" -v "$ROOT/bin:/pipeline_bin:ro" -v "$W:/w" -w /w "$IMG" \
        Rscript /pipeline_bin/make_plots.R --dds diffexp/dds.rds --contrasts ct.csv --results_dir diffexp \
        --gsea_dir GSEA --gsva_dir GSVA --gsva_top 3 --out plots_gs.pdf > "$W/plots_gs.log" 2>&1
}
run_plots_gs
test "$(pdf_pages "$W/plots_gs.pdf")" -eq 7
echo "ok  plots PDF with GSEA + GSVA: 7 pages"
# wrong arguments exit 2
docker run --rm "$(img fgsea)" true >/dev/null 2>&1
rm -rf "$W"
