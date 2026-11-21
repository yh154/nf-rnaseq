// Alignment + quantification (per-sample RSEM files stay in work/: inputs of MERGE_RSEM and Summary.html;
// only the sorted BAMs are published): rsem-calculate-expression --star.
// Lanes are passed as comma-separated lists (no concatenation). The STAR genome BAM is
// unsorted, so it is sorted and indexed here; RSEM's own BAM outputs are switched off.
process STAR_RSEM {
    tag "${sample}"
    label 'star_rsem'
    label 'process_high'

    publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.genome.sorted.bam*'

    input:
    tuple val(sample), val(strandedness), path(r1, stageAs: 'r1/*'), path(r2, stageAs: 'r2/*')
    tuple path(ref_dir), val(prefix)

    output:
    tuple val(sample), path("${sample}.genes.results"), path("${sample}.isoforms.results"), emit: results
    tuple val(sample), path("${sample}.log"), emit: star_log
    path "${sample}.stat", emit: stat
    tuple val(sample), path("${sample}.genome.sorted.bam"), path("${sample}.genome.sorted.bam.bai"), optional: true, emit: bam

    script:
    def r1f   = r1 instanceof List ? r1 : [r1]
    def r2f   = r2 instanceof List ? r2 : (r2 ? [r2] : [])
    def paired = r2f.size() > 0
    def r1l   = r1f.collect { "${it.name}" }.join(',')
    def r2l   = r2f.collect { "${it.name}" }.join(',')
    def gz    = r1f.every { it.name.endsWith('.gz') } ? '--star-gzipped-read-file' : ''
    def pe    = paired ? '--paired-end' : ''
    def bam   = params.no_bam ? '' : '--star-output-genome-bam'
    def sortmem = Math.max(256, (task.memory.toMega() / (2 * task.cpus)) as int)
    """
    trap 'rm -rf rsem_tmp' EXIT
    rsem-calculate-expression --star ${gz} ${pe} \\
        --strandedness ${strandedness} ${bam} --no-bam-output \\
        --temporary-folder rsem_tmp -p ${task.cpus} \\
        ${r1l} ${r2l} ${ref_dir}/${prefix} ${sample} > rsem.stdout 2> rsem.stderr \\
        || { tail -n 40 rsem.stderr >&2; tail -n 15 rsem.stdout >&2; exit 1; }

    if [[ -e ${sample}.STAR.genome.bam ]]; then
        samtools sort -@ ${task.cpus} -m ${sortmem}M -T sorttmp -o ${sample}.genome.sorted.bam ${sample}.STAR.genome.bam
        samtools index -@ ${task.cpus} ${sample}.genome.sorted.bam
        samtools quickcheck -v ${sample}.genome.sorted.bam
        rm -f ${sample}.STAR.genome.bam
    fi
    """
}
