#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * nf-rnaseq: bulk RNA-seq, STAR + RSEM, DESeq2, optional GSEA / GSVA.
 * Phase 1: parameter checks, --help, input validation (later phases add the analysis).
 */

include { VALIDATE_INPUTS } from './modules/validate_inputs.nf'
include { FASTQC } from './modules/fastqc.nf'
include { STAR_RSEM } from './modules/star_rsem.nf'
include { MERGE_RSEM } from './modules/expression.nf'
include { DESEQ2 } from './modules/deseq2.nf'
include { PLOTS } from './modules/plots.nf'
include { VALIDATE_GMT; GSEA; GSVA } from './modules/gene_sets.nf'
include { SUMMARY } from './modules/summary.nf'
include { BUILD_ANNOTATION } from './modules/annotation.nf'
include { DETECT_READ_LENGTH; TOOL_VERSIONS; BUILD_REFERENCE; CHECK_REFERENCE; HASH_INPUTS; VERSIONS } from './modules/reference.nf'

def helpMessage() {
    return """
    Usage: nextflow run . -profile <executor>,<container> --input samples.csv --contrast contrasts.csv \\
               --outdir results (--rsem_ref DIR | --gfasta genome.fa) --gtf genome.gtf [options]

    nf-rnaseq ${workflow.manifest.version}: bulk RNA-seq with STAR + RSEM (rsem-calculate-expression --star),
    DESeq2 and optional GSEA / GSVA.

    Input / output
      --input FILE            Samplesheet CSV: sample,fastq_1,fastq_2,strandedness,condition[,batch]  [required]
                              (empty fastq_2 = single-end; repeated sample rows = lanes;
                              strandedness: forward | reverse | none; absolute paths are used as
                              given, relative paths are resolved against --fastq_dir if set,
                              otherwise against the samplesheet's directory)
      --fastq_dir DIR         Directory of the FASTQ files named (relative) in the samplesheet  [optional]
      --contrast FILE         Desired contrasts  [required]. Either a list of pairs, one per line
                              "treatment,control" (comma or tab; optional third column = variable,
                              default condition; '#' comments), or a CSV id,variable,treatment,control.
                              log2FC = treatment/control
      --outdir DIR            Output directory  [required]

    Reference
      --rsem_ref DIR          Pre-built RSEM+STAR reference directory (skips building)
      --gfasta FILE           Genome FASTA (required when --rsem_ref is not given)
      --gtf FILE              Genome GTF (always required: gene names and biotypes)
      --star_sjdboverhang N   STAR sjdbOverhang when building: 'auto' = read length - 1 [auto]

    DESeq2
      --design STR            Additive design formula ["~ condition"], e.g. "~ batch + condition"
      --padj FLOAT            Adjusted p-value cutoff, also the independentFiltering alpha [0.05]
      --lfc FLOAT             |log2FC| threshold for the volcano plot only [1]
      --min_count INT         Pre-filter: minimum total count over all samples [10]

    GSEA / GSVA / plots
      --gmt FILE[,FILE]       Gene-set file(s) in GMT format (e.g. MSigDB Hallmark/KEGG/GO, or your own;
                              one set per line: name, description, gene1, gene2, ...). Separate several
                              files with commas, no spaces: --gmt h.all.gmt,c2.kegg.gmt
                              Giving this option turns on GSEA (fgsea, ranked by the DESeq2 Wald
                              statistic) and GSVA, run for every contrast x GMT file.
                              Without it neither runs. No gene sets are bundled. Genes must use the
                              species and identifiers of your GTF (see --gmt_id_type). [not set]
      --gmt_id_type TYPE      Identifiers used in the GMT genes: gene_name (symbols, e.g. TP53; MSigDB
                              files) or gene_id (e.g. ENSG..., as in the GTF gene_id)  [gene_name]
      --min_set_size INT      Smallest gene set analysed by GSEA / GSVA, counted after matching to the
                              genes in the data [5]
      --max_set_size INT      Largest gene set analysed by GSEA / GSVA [500]
      --gsva_top INT          Gene sets shown in the GSVA plot [10]
      --gsva_input TYPE       vst | log2tpm  [vst]
      --pca_ntop INT          Top variable genes for PCA [500]
      --run_date DDMMYY       Date in the plot PDF name [today]

    Misc
      --no_bam                Do not write genome-coordinate BAMs
      --queue NAME            Slurm partition (with -profile slurm)
      --help                  Show this message and exit
      --version               Print the version and exit

    Profiles (combine one executor with one container profile)
      local | slurm           executors
      docker | apptainer      containers (images pinned in conf/containers.config)
      test                    synthetic test data (tests/data)

    Examples
      Local Docker (development):
        nextflow run . -profile local,docker --input samples.csv --contrast contrasts.csv \\
            --gfasta genome.fa --gtf genome.gtf --outdir results
      Slurm + Apptainer:
        nextflow run . -profile slurm,apptainer --queue <partition> -c site.config \\
            --input samples.csv --contrast contrasts.csv --gfasta genome.fa --gtf genome.gtf \\
            --outdir results
      Pre-built reference:
        nextflow run . -profile slurm,apptainer --rsem_ref /refs/rsem_star --gtf genome.gtf \\
            --input samples.csv --contrast contrasts.csv --outdir results
    """.stripIndent()
}

def fail(String msg) {
    log.error "${msg}\n\nRun with --help for usage."
    System.exit(1)
}

def checkParams() {
    def known = [
        'input', 'contrast', 'fastq_dir', 'outdir', 'gfasta', 'gtf', 'rsem_ref', 'star_sjdboverhang',
        'design', 'padj', 'lfc', 'min_count', 'gmt', 'gmt_id_type',
        'gsva_top', 'gsva_input', 'min_set_size', 'max_set_size', 'pca_ntop', 'run_date', 'no_bam', 'queue',
        'help', 'version'
    ] as Set
    def unknown = params.keySet().findAll { !(it in known) }
    if (unknown) fail("Unknown parameter(s): ${unknown.collect { '--' + it }.join(', ')}")

    ['input', 'contrast', 'fastq_dir', 'outdir', 'gfasta', 'gtf', 'rsem_ref', 'gmt', 'design', 'queue',
     'run_date', 'gmt_id_type', 'gsva_input', 'star_sjdboverhang'].each {
        if (params[it] instanceof Boolean) fail("--${it} needs a value")
    }

    def missing = ['input', 'outdir', 'gtf'].findAll { !params[it] }
    if (!params.contrast) missing << 'contrast'
    if (!params.rsem_ref && !params.gfasta) missing << 'gfasta (or --rsem_ref)'
    if (missing) fail("Missing required parameter(s): ${missing.collect { '--' + it }.join(', ')}")

    ['input', 'contrast', 'gtf'].each { if (params[it] && !file(params[it]).exists()) fail("--${it}: file not found: ${params[it]}") }
    if (params.fastq_dir && !file(params.fastq_dir).isDirectory()) fail("--fastq_dir: directory not found: ${params.fastq_dir}")
    if (params.gfasta && !file(params.gfasta).exists()) fail("--gfasta: file not found: ${params.gfasta}")
    if (params.rsem_ref && !file(params.rsem_ref).isDirectory()) fail("--rsem_ref: directory not found: ${params.rsem_ref}")
    if (params.gmt) {
        params.gmt.toString().split(',').each { if (!file(it.trim()).exists()) fail("--gmt: file not found: ${it}") }
    }

    def sj = params.star_sjdboverhang.toString()
    if (sj != 'auto' && !(sj ==~ /[1-9][0-9]*/)) fail("--star_sjdboverhang must be 'auto' or a positive integer (got '${sj}')")
    if (!(params.padj instanceof Number) || params.padj <= 0 || params.padj > 1) fail("--padj must be a number in (0, 1]")
    if (!(params.lfc instanceof Number) || params.lfc < 0) fail("--lfc must be a number >= 0")
    if (!(params.min_count instanceof Integer) || params.min_count < 0) fail("--min_count must be an integer >= 0")
    if (!(params.gsva_top instanceof Integer) || params.gsva_top < 1) fail("--gsva_top must be an integer >= 1")
    if (!(params.pca_ntop instanceof Integer) || params.pca_ntop < 2) fail("--pca_ntop must be an integer >= 2")
    if (!(params.min_set_size instanceof Integer) || params.min_set_size < 1) fail("--min_set_size must be an integer >= 1")
    if (!(params.max_set_size instanceof Integer) || params.max_set_size < params.min_set_size) fail("--max_set_size must be an integer >= --min_set_size")
    if (!(params.gmt_id_type in ['gene_name', 'gene_id'])) fail("--gmt_id_type must be gene_name or gene_id")
    if (!(params.gsva_input in ['vst', 'log2tpm'])) fail("--gsva_input must be vst or log2tpm")
    // Nextflow turns a numeric CLI value such as 010101 into the number 10101: restore the zeros
    if (!(params.run_date.toString().padLeft(6, '0') ==~ /\d{6}/)) fail("--run_date must be DDMMYY (got '${params.run_date}')")
    if (!params.design.toString().trim().startsWith('~')) fail("--design must start with '~'")
}

workflow {
    main:
    if (params.version) {
        println "nf-rnaseq ${workflow.manifest.version}"
        System.exit(0)
    }
    if (params.help) {
        println helpMessage()
        System.exit(0)
    }
    checkParams()

    def samplesheet = file(params.input)
    VALIDATE_INPUTS(
        samplesheet,
        file(params.contrast),
        params.design,
        (params.fastq_dir ? file(params.fastq_dir) : samplesheet.parent).toAbsolutePath().toString()
    )

    // FASTQ files are checked on the host (they may not be mounted in the container).
    VALIDATE_INPUTS.out.samplesheet
        .splitCsv(header: true)
        .map { row ->
            [row.fastq_1, row.fastq_2].findAll { it }.each { f ->
                if (!file(f).exists()) fail("FASTQ not found: ${f}")
            }
            row
        }
        .set { ch_rows }

    // ---- reference and annotation
    BUILD_ANNOTATION(file(params.gtf))
    TOOL_VERSIONS()

    // read length of the first sample (all FASTQ mates): sjdbOverhang = max read length - 1
    VALIDATE_INPUTS.out.samplesheet
        .splitCsv(header: true)
        .first()
        .map { row -> [row.fastq_1, row.fastq_2].findAll { it }.collect { file(it) } }
        .set { ch_first_fastqs }
    DETECT_READ_LENGTH(ch_first_fastqs)

    if (params.rsem_ref) {
        ch_ref_dir = Channel.value(file(params.rsem_ref))
        ch_hash = Channel.value([[file(params.gtf), file(params.rsem_ref)], [['gtf', file(params.gtf).name], ['rsem_ref', file(params.rsem_ref).name]]])
    } else {
        ch_overhang = params.star_sjdboverhang.toString() == 'auto'
            ? DETECT_READ_LENGTH.out.length.map { (it.text.trim() as int) - 1 }
            : Channel.value(params.star_sjdboverhang as int)
        BUILD_REFERENCE(file(params.gfasta), file(params.gtf), ch_overhang)
        ch_ref_dir = BUILD_REFERENCE.out.dir
        ch_hash = Channel.value([[file(params.gfasta), file(params.gtf)], [['gfasta', file(params.gfasta).name], ['gtf', file(params.gtf).name]]])
    }
    CHECK_REFERENCE(ch_ref_dir, BUILD_ANNOTATION.out.annotation, BUILD_ANNOTATION.out.contigs,
                    TOOL_VERSIONS.out.star_version, DETECT_READ_LENGTH.out.length)
    HASH_INPUTS(ch_hash.map { it[0] }, ch_hash.map { it[1] })

    // [reference dir, RSEM prefix] for the alignment step (Phase 3)
    CHECK_REFERENCE.out.prefix
        .combine(ch_ref_dir)
        .map { pf, dir -> [dir, pf.text.trim()] }
        .first()
        .set { ch_ref }

    // ---- per-FASTQ QC
    ch_rows
        .flatMap { r -> [[r.sample, file(r.fastq_1)]] + (r.fastq_2 ? [[r.sample, file(r.fastq_2)]] : []) }
        .set { ch_fastq }
    FASTQC(ch_fastq)

    // ---- per-sample alignment + quantification (lanes grouped, order kept)
    ch_rows
        .map { r -> [r.sample, r.strandedness, file(r.fastq_1), r.fastq_2 ? file(r.fastq_2) : null] }
        .groupTuple()
        .map { sample, strand, r1, r2 -> [sample, strand[0], r1, r2.findAll { it }] }
        .set { ch_samples }
    STAR_RSEM(ch_samples, ch_ref)

    // ---- expression tables and differential expression
    STAR_RSEM.out.results.map { it[1] }.collect().set { ch_genes_results }
    MERGE_RSEM(VALIDATE_INPUTS.out.samplesheet, BUILD_ANNOTATION.out.annotation, ch_genes_results)
    DESEQ2(MERGE_RSEM.out.counts, VALIDATE_INPUTS.out.samplesheet, VALIDATE_INPUTS.out.contrasts,
           params.design, params.padj, params.min_count)

    // ---- plots (one PDF) and Summary.html
    DESEQ2.out.results.collect().ifEmpty([]).set { ch_deseq_results }

    // ---- optional gene-set analyses (--gmt): GSEA on the Wald statistic, GSVA on VST / log2 TPM
    if (params.gmt) {
        ch_gmt = Channel.value(params.gmt.toString().split(',').collect { file(it.trim()) })
        VALIDATE_GMT(ch_gmt, BUILD_ANNOTATION.out.annotation, params.gmt_id_type, params.min_set_size)
        GSEA(ch_deseq_results, VALIDATE_INPUTS.out.contrasts, ch_gmt, VALIDATE_GMT.out.report,
             params.gmt_id_type, params.min_set_size, params.max_set_size)
        ch_gsva_expr = params.gsva_input == 'vst' ? DESEQ2.out.vst : MERGE_RSEM.out.tpm
        GSVA(ch_gsva_expr, VALIDATE_INPUTS.out.samplesheet, VALIDATE_INPUTS.out.contrasts, ch_gmt,
             VALIDATE_GMT.out.report, params.gsva_input, params.gmt_id_type, params.min_set_size, params.max_set_size)
        VALIDATE_GMT.out.report.map { it.text }.filter { it.contains('WARNING') }
            .subscribe { txt -> log.warn "Gene sets: " + txt.readLines().findAll { it.startsWith('WARNING') }.join('; ') }
        ch_gsea_tables = GSEA.out.tables.collect().ifEmpty([])
        ch_gsva_tables = GSVA.out.tables.collect().ifEmpty([])
    } else {
        ch_gsea_tables = Channel.value([])
        ch_gsva_tables = Channel.value([])
    }
    // ---- pipeline_info/versions.yml
    def images = file("${projectDir}/conf/containers.config").readLines()
        .findAll { it ==~ /\s+\/\/ \w+: .*/ }.collect { "  ${it.trim().substring(3)}" }
    def vHeader = (["pipeline: nf-rnaseq ${workflow.manifest.version}", "nextflow: ${nextflow.version}", 'containers:'] +
        images + ['tools:']).join('\n').replace("'", '')
    ch_versions = TOOL_VERSIONS.out.versions.mix(DESEQ2.out.versions)
    if (params.gmt) ch_versions = ch_versions.mix(GSEA.out.versions, GSVA.out.versions)
    VERSIONS(ch_versions.collect(), vHeader)

    PLOTS(DESEQ2.out.dds, VALIDATE_INPUTS.out.contrasts, ch_deseq_results, ch_gsea_tables, ch_gsva_tables,
          params.padj, params.lfc, params.pca_ntop, params.gsva_top, params.gmt ? true : false)

    STAR_RSEM.out.results.join(STAR_RSEM.out.star_log)
        .map { sample, genes, iso, log -> [genes, log] }
        .flatten().collect().set { ch_rsem_files }
    def footer = "Generated by nf-rnaseq ${workflow.manifest.version} on ${new Date().format('yyyy-MM-dd')}".toString().replace("'", '')
    def info = [
        "Pipeline: nf-rnaseq ${workflow.manifest.version}",
        "Design: ${params.design}",
        "Reference: ${params.rsem_ref ? 'pre-built ' + file(params.rsem_ref).name : 'built from ' + file(params.gfasta).name}",
        "Annotation: ${file(params.gtf).name}",
        "padj / min_count: ${params.padj} / ${params.min_count}"
    ].join('\n').replace("'", '')
    VALIDATE_INPUTS.out.samplesheet.splitCsv(header: true).first().map { it.fastq_2 ? false : true }.set { ch_single_end }
    SUMMARY(FASTQC.out.zip.collect(), ch_rsem_files, BUILD_ANNOTATION.out.annotation,
            VALIDATE_INPUTS.out.samplesheet, VALIDATE_INPUTS.out.warnings, ch_single_end, footer, info)
    SUMMARY.out.assigned_warnings.map { it.text }.filter { it.trim() }
        .subscribe { txt -> log.warn "Few reads assigned to genes:\n" + txt.trim() }

    onComplete:
    if (!params.help && !params.version && params.outdir) {
        def infoDir = file("${params.outdir}/pipeline_info")
        infoDir.mkdirs()
        def p = new TreeMap(params.findAll { k, v -> !k.contains('-') && !(k ==~ /.*[A-Z].*/) })
        file("${infoDir}/params.json").text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(
            [pipeline: workflow.manifest.name, version: workflow.manifest.version,
             nextflow: nextflow.version.toString(), success: workflow.success, params: p]))
    }
}
