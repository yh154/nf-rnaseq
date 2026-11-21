// One multi-page PDF. The date (--run_date) is used ONLY in the published file name, never
// as a process input, so a later -resume does not redraw the plots.
process PLOTS {
    label 'deseq2'
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> "plot/plots_${params.run_date.toString().padLeft(6, '0')}.pdf" }

    input:
    path dds
    path contrasts
    path results, stageAs: 'results/*'
    path gsea, stageAs: 'gsea/*'
    path gsva, stageAs: 'gsva/*'
    val padj
    val lfc
    val pca_ntop
    val gsva_top
    val has_gmt

    output:
    path 'plots.pdf', emit: pdf

    script:
    def gs = has_gmt ? "--gsea_dir gsea --gsva_dir gsva --gsva_top ${gsva_top}" : ''
    """
    mkdir -p results
    make_plots.R --dds ${dds} --contrasts ${contrasts} --results_dir results \\
        --padj ${padj} --lfc ${lfc} --pca_ntop ${pca_ntop} ${gs} --out plots.pdf
    """
}
