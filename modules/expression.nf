// Gene-level expression tables from the per-sample RSEM results.
process MERGE_RSEM {
    label 'python'
    label 'process_low'

    publishDir "${params.outdir}/expression", mode: 'copy', saveAs: { fn -> fn.tokenize('/')[-1] }

    input:
    path samplesheet
    path annotation
    path genes_results

    output:
    path 'expr/RSEM_expected_counts.txt', emit: counts
    path 'expr/TPM.txt',                  emit: tpm
    path 'expr/FPKM.txt',                 emit: fpkm

    script:
    """
    merge_rsem.py --samplesheet ${samplesheet} --annotation ${annotation} --outdir expr ${genes_results}
    """
}
