// One DESeq2 fit (vst_counts.txt stays in work/: GSVA input only, not published) on round(expected_count); one results() per contrast.
process DESEQ2 {
    label 'deseq2'
    label 'process_medium'

    publishDir "${params.outdir}", mode: 'copy', saveAs: { fn ->
        def n = fn.tokenize('/')[-1]
        n == 'DESeq2_normalized_counts.txt' ? "expression/${n}" :
        n in ['versions_deseq2.yml', 'vst_counts.txt'] ? null : "diffexp/${n}"
    }

    input:
    path counts
    path samplesheet
    path contrasts
    val design
    val padj
    val min_count

    output:
    path 'diffexp/dds.rds',                       emit: dds
    path 'diffexp/Deseq2_res_*.csv',              emit: results, optional: true
    path 'diffexp/DESeq2_normalized_counts.txt',  emit: normalized
    path 'diffexp/vst_counts.txt',                emit: vst
    path 'diffexp/versions_deseq2.yml',           emit: versions

    script:
    """
    run_deseq2.R --counts ${counts} --samplesheet ${samplesheet} --contrasts ${contrasts} \\
        --design '${design}' --padj ${padj} --min_count ${min_count} --outdir diffexp
    """
}
