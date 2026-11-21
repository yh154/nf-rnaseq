process VALIDATE_INPUTS {
    label 'python'
    label 'process_low'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { fn -> fn in ['validation_warnings.txt'] ? fn : null }

    input:
    path samplesheet
    path contrasts
    val design
    val base_dir

    output:
    path 'samplesheet.valid.csv', emit: samplesheet
    path 'contrasts.valid.csv',   emit: contrasts
    path 'validation_warnings.txt', emit: warnings

    script:
    """
    validate_samplesheet.py \\
        --samplesheet ${samplesheet} --contrasts ${contrasts} \\
        --design '${design}' --base_dir '${base_dir}' --no_file_check \\
        --out_samplesheet samplesheet.valid.csv \\
        --out_contrasts contrasts.valid.csv \\
        --out_warnings validation_warnings.txt
    """
}
