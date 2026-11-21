// Summary.html from the FastQC reports and the STAR / RSEM results (bash + awk + unzip).
process SUMMARY {
    label 'python'
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'Summary.html'

    input:
    path fastqc_zips, stageAs: 'fastqc/*'
    path rsem_files,  stageAs: 'rsem/*'
    path annotation
    path samplesheet
    path warnings
    val single_end
    val footer
    val info

    output:
    path 'Summary.html', emit: html
    path 'assigned_warnings.txt', emit: assigned_warnings

    script:
    def se = single_end ? '--single-end' : ''
    """
    awk -F, 'NR > 1 && !seen[\$1]++ { print \$1 }' ${samplesheet} > sample_order.txt
    printf '%s\\n' '${info}' > run_info.txt
    check_assigned.py --rsem_dir rsem --out assigned_warnings.txt
    cat ${warnings} assigned_warnings.txt > all_warnings.txt
    summaryHTML.sh ${se} --annotation ${annotation} --sample_order sample_order.txt \\
        --info run_info.txt --warnings all_warnings.txt --footer '${footer}' fastqc rsem .
    """
}
