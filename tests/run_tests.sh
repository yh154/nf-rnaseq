#!/usr/bin/env bash
# Single test entry point.
#   tests/run_tests.sh [--no-integration] [--profile PROFILES]
# Steps: shellcheck (if installed), synthetic data, python unit tests, --help tests,
# integration (nextflow, default profile "test,local"; use --profile test,docker once the
# images are built). Prints a PASS/FAIL/SKIP line per step; exit 1 if any step failed.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 1
INTEGRATION=1
PROFILE="test,local"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-integration) INTEGRATION=0 ;;
        --profile) PROFILE="$2"; shift ;;
        -h|--help)
            sed -n '2,7p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            echo "Usage: tests/run_tests.sh [--no-integration] [--profile PROFILES]"
            exit 0 ;;
        *) echo "Usage: tests/run_tests.sh [--no-integration] [--profile PROFILES]" >&2; exit 2 ;;
    esac
    shift
done

FAILED=0
step() { # name, command...
    local name="$1"; shift
    if "$@" >"$WORK/$name.log" 2>&1; then echo "PASS  $name"
    else echo "FAIL  $name  (log: $WORK/$name.log)"; tail -15 "$WORK/$name.log" | sed 's/^/      /'; FAILED=1; fi
}
skip() { echo "SKIP  $1 ($2)"; }
WORK="$(mktemp -d "${TMPDIR:-/tmp}/nfrnaseq_tests.XXXXXX")"

# Lint with a local binary, else the pinned image (development tool only, not a pipeline container)
SHELLCHECK_IMAGE="koalaman/shellcheck@sha256:2097951f02e735b613f4a34de20c40f937a6c8f18ecb170612c88c34517221fb"
if command -v shellcheck >/dev/null; then
    step shellcheck shellcheck bin/*.sh tests/*.sh
elif command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
    step shellcheck docker run --rm -v "$ROOT:/mnt" -w /mnt "$SHELLCHECK_IMAGE" bin/*.sh tests/*.sh
else
    skip shellcheck "neither shellcheck nor docker available - untested"
fi

step synthetic_data python3 -I tests/make_synthetic.py
step unit_validate python3 -W ignore::ResourceWarning -I tests/test_validate.py
step unit_phase2 python3 -W ignore::ResourceWarning -I tests/test_phase2.py
step unit_merge python3 -W ignore::ResourceWarning -I tests/test_merge.py
step unit_summary python3 -W ignore::ResourceWarning -I tests/test_summary.py
step unit_gmt python3 -W ignore::ResourceWarning -I tests/test_gmt.py
step unit_assigned python3 -W ignore::ResourceWarning -I tests/test_assigned.py
step help_scripts bash tests/test_help.sh
if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
    step phase4_fast bash tests/test_phase4_fast.sh
else
    skip phase4_fast "docker not available"
fi
step pin_script_help python3 -I tools/pin_containers.py --help

if [[ $INTEGRATION -eq 1 ]]; then
    if command -v nextflow >/dev/null; then
        NF=nextflow
    else
        NF=""
    fi
    if [[ -z "$NF" ]]; then
        skip integration "nextflow not on PATH"
    else
        step profiles_resolve bash tests/test_profiles.sh "$NF"
        step main_help bash -c "$NF run . --help | grep -q Usage"
        step main_missing_params bash -c "! $NF run . >/dev/null 2>&1"
        step main_unknown_param bash -c "! $NF run . -profile $PROFILE --outDir x >/dev/null 2>&1"
        NFR="$NF run . -profile $PROFILE -ansi-log false"
        # shellcheck disable=SC2086  # NFR is a command line, split on purpose
        step validate_and_build $NFR --outdir "$WORK/out"
        R="$WORK/out/reference"
        step build_outputs bash -c "
            test -d $R/rsem_star && test -s $R/annotation.tsv &&
            test \$(wc -l < $R/annotation.tsv) -eq 41 &&
            head -1 $R/annotation.tsv | grep -qx 'gene_id.gene_name.gene_biotype' &&
            grep -q '^sjdbOverhang.74\$' $R/rsem_star/genomeParameters.txt &&
            grep -q sha256 $WORK/out/pipeline_info/inputs_sha256.json &&
            test -s $WORK/out/pipeline_info/params.json && test -s $WORK/out/pipeline_info/trace.txt"
        step versions_yml bash -c "
            grep -q '^pipeline: nf-rnaseq' $WORK/out/pipeline_info/versions.yml &&
            grep -q '  STAR: 2.7.10a' $WORK/out/pipeline_info/versions.yml &&
            test \"\$(ls $WORK/out/pipeline_info | paste -sd, -)\" = 'inputs_sha256.json,params.json,trace.txt,validation_warnings.txt,versions.yml' &&
            grep -q '  DESeq2: ' $WORK/out/pipeline_info/versions.yml &&
            grep -q '  fgsea: ' $WORK/out/pipeline_info/versions.yml &&
            grep -q '  GSVA: ' $WORK/out/pipeline_info/versions.yml &&
            grep -q '  star_rsem: mulled-v2' $WORK/out/pipeline_info/versions.yml"
        step resume_reuses_reference bash -c "$NFR --outdir $WORK/out -resume > $WORK/resume.out 2>&1; awk -F'\t' '\$4==\"BUILD_REFERENCE\" && \$5==\"CACHED\"' $WORK/out/pipeline_info/trace.txt | grep -q ."
        step phase3_pe python3 -I tests/check_phase3.py "$WORK/out" tests/data/samples_pe.csv
        DESEQ_IMG="$(grep 'withLabel: deseq2 ' conf/containers.config | sed "s/.*container = '\(.*\)'.*/\1/")"
        step phase4_tables python3 -I tests/check_phase4.py tables "$WORK/out/expression" "$WORK/out" tests/data/samples_pe.csv
        step phase4_deseq python3 -I tests/check_phase4.py deseq "$WORK/out/diffexp" "$WORK/out/expression" tests/data/samples_pe.csv
        step phase4_golden python3 -I tests/check_phase4.py golden "$WORK/out/expression"
        step phase4_dds_rds docker run --rm -v "$ROOT/tests:/t:ro" -v "$WORK/out/diffexp:/d:ro" "$DESEQ_IMG" \
            Rscript /t/check_dds.R /d/dds.rds "~ batch + condition" 6
        step phase5_summary_and_pdf python3 -I tests/check_phase5.py "$WORK/out" tests/data/samples_pe.csv 1 "$(date +%d%m%y)" 2
        # --run_date only renames the published PDF: -resume must not redraw the plots
        step run_date_does_not_rerun_plots bash -c "
            $NFR --outdir $WORK/out_rd --run_date 010101 -resume > $WORK/rd.out 2>&1;
            awk -F'\t' '\$4==\"PLOTS\" && \$5==\"CACHED\"' $WORK/out_rd/pipeline_info/trace.txt | grep -q . &&
            test -s $WORK/out_rd/plot/plots_010101.pdf"
        step phase6_gsea_gsva python3 -I tests/check_phase6.py "$WORK/out/gsea" "$WORK/out/gsva"
        step gmt_wrong_id_type_fails bash -c "
            ! $NFR --outdir $WORK/out_gid --rsem_ref $R/rsem_star --gmt_id_type gene_id > $WORK/gid.out 2>&1 &&
            grep -q 'wrong --gmt_id_type' $WORK/gid.out"
        step no_leftover_tmp bash -c "! find work -type d \\( -name rsem_tmp -o -name '*_STARtmp' \\) | grep -q ."
        step prebuilt_reference bash -c "
            $NFR --outdir $WORK/out_pre --rsem_ref $R/rsem_star > $WORK/pre.out 2>&1;
            ! grep -q 'BUILD_REFERENCE' $WORK/pre.out && grep -q 'CHECK_REFERENCE' $WORK/pre.out &&
            grep -q rsem_ref $WORK/out_pre/pipeline_info/inputs_sha256.json"
        step prebuilt_same_counts_as_built python3 -I tests/compare_results.py "$WORK/out" "$WORK/out_pre"
        # --contrast: the same contrast given as a plain pair list (overrides the full-format contrast file of the test profile)
        step phase3_se_prebuilt_contrast_pairs bash -c "
            $NFR --outdir $WORK/out_se --rsem_ref $R/rsem_star --input tests/data/samples_se.csv \\
                --contrast tests/data/contrast_pairs.txt > $WORK/se.out 2>&1 &&
            test -s $WORK/out_se/diffexp/Deseq2_res_treated_vs_control.csv &&
            python3 -I tests/check_phase3.py $WORK/out_se tests/data/samples_se.csv"
        step contrast_file_not_found bash -c "
            ! $NFR --outdir $WORK/out_cnf --contrast /no/such/pairs.txt > $WORK/cnf.out 2>&1 &&
            grep -q -- '--contrast: file not found' $WORK/cnf.out"
        step fastq_dir_not_found bash -c "
            ! $NFR --outdir $WORK/out_fd --fastq_dir /no/such/dir > $WORK/fd.out 2>&1 &&
            grep -q -- '--fastq_dir: directory not found' $WORK/fd.out"
        step contrast_required bash -c "
            ! $NF run . -profile local --input x --outdir y --gtf z --gfasta w > $WORK/creq.out 2>&1 &&
            grep -q -- '--contrast' $WORK/creq.out"
        # four samples (2 vs 2), no BAM, overridden sjdbOverhang (builds its own reference)
        grep -E '^(sample|control_rep1|control_rep2|treated_rep1|treated_rep2),' tests/data/samples_pe.csv | sed 's#,fastq/#,'"$ROOT"'/tests/data/fastq/#g' > "$WORK/four.csv"
        step no_bam_and_overhang_override bash -c "
            $NFR --outdir $WORK/out_nb --input $WORK/four.csv --design '~ condition' \\
                --no_bam --star_sjdboverhang 49 --gsva_input log2tpm > $WORK/nb.out 2>&1 &&
            test -s $WORK/out_nb/gsva/GSVA_scores_test.csv &&
            test -s $WORK/out_nb/diffexp/Deseq2_res_treated_vs_control.csv &&
            test -z \"\$(ls $WORK/out_nb/bam 2>/dev/null)\" &&
            grep -q '^sjdbOverhang.49\$' $WORK/out_nb/reference/rsem_star/genomeParameters.txt &&
            python3 -I tests/check_phase3.py $WORK/out_nb $WORK/four.csv --no-bam"
        # an unreplicated sample (only one with condition=other) is skipped from the entire analysis:
        # not aligned, not in any table, and the other samples' results are unchanged (same size factors)
        mkdir -p "$WORK/lone"
        ln -sf "$ROOT/tests/data/fastq/treated_rep1_pe_R1.fq.gz" "$WORK/lone/lone_R1.fq.gz"
        ln -sf "$ROOT/tests/data/fastq/treated_rep1_pe_R2.fq.gz" "$WORK/lone/lone_R2.fq.gz"
        { sed 's#,fastq/#,'"$ROOT"'/tests/data/fastq/#g' tests/data/samples_pe.csv;
          echo "lone,$WORK/lone/lone_R1.fq.gz,$WORK/lone/lone_R2.fq.gz,reverse,other,b1"; } > "$WORK/with_lone.csv"
        step unreplicated_sample_skipped bash -c "
            $NFR --outdir $WORK/out_lone --rsem_ref $R/rsem_star --input $WORK/with_lone.csv > $WORK/lone.out 2>&1 &&
            grep -q 'sample \"lone\" skipped' $WORK/out_lone/pipeline_info/validation_warnings.txt &&
            ! grep -qw lone $WORK/out_lone/expression/RSEM_expected_counts.txt &&
            ! grep -qw lone $WORK/out_lone/expression/DESeq2_normalized_counts.txt &&
            ! grep -q 'STAR_RSEM (lone)' $WORK/out_lone/pipeline_info/trace.txt &&
            ! grep -qw lone $WORK/out_lone/Summary.html &&
            cmp $WORK/out/expression/DESeq2_normalized_counts.txt $WORK/out_lone/expression/DESeq2_normalized_counts.txt &&
            cmp $WORK/out/diffexp/Deseq2_res_treated_vs_control.csv $WORK/out_lone/diffexp/Deseq2_res_treated_vs_control.csv"
        # strandedness declared as the opposite of the library: RSEM must fail with a clear message
        sed 's/,reverse,control,b1/,forward,control,b1/' "$WORK/four.csv" > "$WORK/wrong_strand.csv"
        step wrong_strandedness_fails bash -c "
            ! $NFR --outdir $WORK/out_ws --rsem_ref $R/rsem_star --input $WORK/wrong_strand.csv --design '~ condition' \\
                > $WORK/ws.out 2>&1 &&
            grep -q 'No valid read to estimate the length distribution' $WORK/ws.out"
        sed 's/^chr/NC_/' tests/data/genes.gtf > "$WORK/renamed.gtf"
        step contig_mismatch_fails bash -c "
            ! $NFR --outdir $WORK/out_bad --rsem_ref $R/rsem_star --gtf $WORK/renamed.gtf > $WORK/bad.out 2>&1 &&
            grep -q 'of GTF contigs are in the reference' $WORK/bad.out"
    fi
else
    skip integration "--no-integration"
fi

if [[ $FAILED -eq 0 ]]; then
    rm -rf "$WORK"      # a passing run leaves several GB of outputs; failures keep them (logs above)
    echo "ALL TESTS PASSED"
else
    echo "SOME TESTS FAILED (outputs kept in $WORK)"
fi
exit $FAILED
