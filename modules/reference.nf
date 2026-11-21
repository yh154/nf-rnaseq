process DETECT_READ_LENGTH {
    label 'python'
    label 'process_low'

    input:
    path fastqs

    output:
    path 'read_length.txt', emit: length

    script:
    """
    max_read_length.py ${fastqs} > read_length.txt
    """
}

process TOOL_VERSIONS {
    label 'star_rsem'
    label 'process_low'

    output:
    path 'star_version.txt',        emit: star_version
    path 'versions_star_rsem.yml',  emit: versions

    script:
    """
    STAR --version > star_version.txt
    {
        echo "STAR: \$(cat star_version.txt)"
        echo "RSEM: \$(rsem-calculate-expression --version | sed 's/Current version: //')"
        echo "samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')"
    } > versions_star_rsem.yml
    """
}

process BUILD_REFERENCE {
    label 'star_rsem'
    label 'process_high'

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta
    path gtf
    val overhang

    output:
    path 'rsem_star', emit: dir

    script:
    """
    FA=${fasta}; GT=${gtf}
    if [[ "\$FA" == *.gz ]]; then gzip -dc "\$FA" > genome.fa; FA=genome.fa; fi
    if [[ "\$GT" == *.gz ]]; then gzip -dc "\$GT" > genome.gtf; GT=genome.gtf; fi
    mkdir rsem_star
    echo "STAR sjdbOverhang = ${overhang}"
    rsem-prepare-reference --gtf "\$GT" --star --star-sjdboverhang ${overhang} \\
        -p ${task.cpus} "\$FA" rsem_star/ref > build_reference.log 2>&1 \\
        || { tail -n 30 build_reference.log >&2; exit 1; }
    """
}

process CHECK_REFERENCE {
    label 'python'
    label 'process_low'

    input:
    path ref_dir
    path annotation
    path gtf_contigs
    path star_version
    path read_length

    output:
    path 'prefix.txt', emit: prefix

    script:
    """
    check_reference.py --ref_dir ${ref_dir} --annotation ${annotation} \\
        --gtf_contigs ${gtf_contigs} --star_version "\$(cat ${star_version})" \\
        --read_length "\$(cat ${read_length})" --out_prefix prefix.txt
    """
}

process HASH_INPUTS {
    label 'python'
    label 'process_low'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path files
    val labels

    output:
    path 'inputs_sha256.json'

    script:
    def pairs = labels.collect { "${it[0]}=${it[1]}" }.join(' ')
    """
    hash_inputs.py --out inputs_sha256.json ${pairs}
    """
}

// pipeline_info/versions.yml: pipeline, Nextflow, container images and the tool versions the processes recorded
process VERSIONS {
    label 'python'
    label 'process_low'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path tool_files
    val header

    output:
    path 'versions.yml'

    script:
    """
    printf '%s\\n' '${header}' > versions.yml
    for f in \$(ls versions_*.yml | sort); do
        sed -e 's/^ *//' -e '/^\$/d' -e 's/^/  /' \$f >> versions.yml
    done
    """
}
