// Gene-set analyses (only with --gmt): GMT validation, GSEA (fgsea) and GSVA.
process VALIDATE_GMT {
    label 'python'
    label 'process_low'

    input:
    path gmts, stageAs: 'gmt/*'
    path annotation
    val id_type
    val min_size

    output:
    path 'gmt_report.txt', emit: report

    script:
    """
    validate_gmt.py --gmt "\$(ls gmt/* | paste -sd, -)" --annotation ${annotation} \\
        --gmt_id_type ${id_type} --min_set_size ${min_size} --out gmt_report.txt
    """
}

process GSEA {
    label 'fgsea'
    label 'process_medium'

    publishDir "${params.outdir}/gsea", mode: 'copy', pattern: 'GSEA/*.csv', saveAs: { fn -> fn.tokenize('/')[-1] }

    input:
    path results, stageAs: 'results/*'
    path contrasts
    path gmts, stageAs: 'gmt/*'
    path gmt_report
    val id_type
    val min_size
    val max_size

    output:
    path 'GSEA/GSEA_*.csv', emit: tables, optional: true
    path 'GSEA/versions_fgsea.yml', emit: versions

    script:
    """
    run_gsea.R --results_dir results --contrasts ${contrasts} --gmt "\$(ls gmt/* | paste -sd, -)" \\
        --gmt_id_type ${id_type} --min_set_size ${min_size} --max_set_size ${max_size} --outdir GSEA
    """
}

process GSVA {
    label 'gsva'
    label 'process_medium'

    publishDir "${params.outdir}/gsva", mode: 'copy', pattern: 'GSVA/*.csv', saveAs: { fn -> fn.tokenize('/')[-1] }

    input:
    path expr
    path samplesheet
    path contrasts
    path gmts, stageAs: 'gmt/*'
    path gmt_report
    val input_type
    val id_type
    val min_size
    val max_size

    output:
    path 'GSVA/*.csv', emit: tables
    path 'GSVA/versions_gsva.yml', emit: versions

    script:
    """
    run_gsva.R --expr ${expr} --input_type ${input_type} --samplesheet ${samplesheet} --contrasts ${contrasts} \\
        --gmt "\$(ls gmt/* | paste -sd, -)" --gmt_id_type ${id_type} \\
        --min_set_size ${min_size} --max_set_size ${max_size} --outdir GSVA
    """
}
