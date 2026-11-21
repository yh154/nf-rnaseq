process BUILD_ANNOTATION {
    label 'python'
    label 'process_low'

    publishDir "${params.outdir}/reference", mode: 'copy', pattern: 'annotation.tsv'

    input:
    path gtf

    output:
    path 'annotation.tsv', emit: annotation
    path 'gtf_contigs.txt', emit: contigs

    script:
    """
    build_annotation.py --gtf ${gtf} --out annotation.tsv --contigs_out gtf_contigs.txt
    """
}
