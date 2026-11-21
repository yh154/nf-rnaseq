process FASTQC {
    tag "${sample}:${fastq.name}"
    label 'fastqc'
    label 'process_low'

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    path '*_fastqc.zip',  emit: zip
    path '*_fastqc.html', emit: html

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${fastq}
    """
}
