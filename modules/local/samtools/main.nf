/*
========================================================================================
    MODULE: SAMTOOLS_SORT
========================================================================================
*/

process SAMTOOLS_SORT {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/samtools:latest'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam
    path  "versions.yml",                           emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // -m: memory per thread (keep below process memory / threads)
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -m 2G \\
        -o ${prefix}.sorted.bam \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: SAMTOOLS_INDEX
========================================================================================
*/

process SAMTOOLS_INDEX {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/samtools:latest'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam}.bai"), emit: bai
    path  "versions.yml",               emit: versions

    script:
    """
    samtools index -@ ${task.cpus} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: SAMTOOLS_MARKDUP
    Notes:
      - samtools markdup operates on name-sorted OR coordinate-sorted BAMs.
        For coordinate-sorted input, use -S (coordinate-sorted mode).
      - --json: write duplicate stats in JSON format (parsed by MultiQC).
      - For ONT WGS from a single flow cell, expected duplicate rate is <1%.
        High duplicate rates (>5%) may indicate PCR amplification artefacts
        or a low-input library prep.
      - Optical duplicates: not applicable to ONT (sequencing is not optical);
        use -d 0 to disable optical distance filtering.
========================================================================================
*/

process SAMTOOLS_MARKDUP {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/samtools:latest'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.markdup.bam"), emit: bam
    tuple val(meta), path("${meta.id}.markdup.json"), emit: stats
    path  "versions.yml",                             emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools markdup \\
        -@ ${task.cpus} \\
        --json \\
        -S \\
        -d 0 \\
        ${bam} \\
        ${prefix}.markdup.bam \\
        2> ${prefix}.markdup.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.markdup.bam
    echo '{}' > ${meta.id}.markdup.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: SAMTOOLS_FLAGSTAT
========================================================================================
*/

process SAMTOOLS_FLAGSTAT {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/flagstat/${meta.id}", mode: 'copy'

    container   'sktrinh12/samtools:latest'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.flagstat"), emit: stats
    path  "versions.yml",                         emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools flagstat \\
        -@ ${task.cpus} \\
        ${bam} \\
        > ${prefix}.flagstat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    echo '{}' > ${meta.id}.flagstat
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: SAMTOOLS_BAM2FASTQ
    Notes:
      - Converts Dorado uBAM to FASTQ while stripping MM/ML tags (FASTQ cannot
        carry modification tags). The uBAM is kept as primary output; FASTQ is
        only produced for tools that require it (e.g., vg giraffe).
      - -T '*' passes all auxiliary tags as comment fields (FASTQ comment format:
        @readname  TAG:TYPE:VALUE) — modkit can parse these but most FASTQ tools
        cannot, so we use -T '' for clean FASTQ output.
      - -0: write all reads to a single interleaved output (not split by flag).
========================================================================================
*/

process SAMTOOLS_BAM2FASTQ {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/samtools:latest'

    input:
    tuple val(meta), path(ubam)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads
    path  "versions.yml",                         emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // -0: single output (not split PE)
    // -T '': no auxiliary tags in FASTQ (clean output for vg giraffe)
    // Pipe through pigz for parallel compression
    """
    samtools bam2fq \\
        -@ ${task.cpus} \\
        -0 /dev/stdout \\
        -T '' \\
        ${ubam} \\
        | gzip -c > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: "stub"
    END_VERSIONS
    """
}
