/*
========================================================================================
    MODULE: NANOPLOT_BAM
    Tool   : NanoPlot (>=1.41.0)
    Purpose: Generate plots and statistics for aligned BAM files.
             Provides N50, mean Q, read length distributions post-alignment.
========================================================================================
*/

process NANOPLOT_BAM {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/nanoplot_bam/${meta.id}", mode: 'copy'

    container   'sktrinh12/nanoplot:latest'

    input:
    tuple val(meta), path(bam)
    val  extra_opts

    output:
    tuple val(meta), path("${meta.id}_nanoplot_bam"), emit: report
    path  "versions.yml",                              emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def opts = extra_opts ?: ''
    """
    NanoPlot \\
        --bam ${bam} \\
        --outdir ${prefix}_nanoplot_bam \\
        --title ${prefix}_bam_qc \\
        --format png \\
        --threads ${task.cpus} \\
        ${opts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_nanoplot_bam
    touch ${prefix}_nanoplot_bam/summary.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: NANOPLOT_RAW
    Tool   : NanoPlot (>=1.41.0)
    Purpose: Generate plots and statistics for pre-basecalled FASTQ input.
========================================================================================
*/

process NANOPLOT_RAW {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/nanoplot_raw/${meta.id}", mode: 'copy'

    container   'sktrinh12/nanoplot:latest'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_nanoplot_raw"), emit: report
    path  "versions.yml",                              emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    NanoPlot \\
        --fastq ${reads} \\
        --outdir ${prefix}_nanoplot_raw \\
        --title ${prefix}_raw_qc \\
        --format png \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_nanoplot_raw
    touch ${prefix}_nanoplot_raw/summary.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: "stub"
    END_VERSIONS
    """
}
