/*
========================================================================================
    MODULE: MOSDEPTH
    Tool   : mosdepth (>=0.3.6)
    Purpose: Fast, accurate genome-wide coverage depth calculation.
             Used as a QC gate: samples below params.min_coverage are flagged.

    Output files:
      *.mosdepth.summary.txt   — per-chromosome + total mean/min/max coverage
      *.mosdepth.global.dist.txt — cumulative coverage distribution
      *.per-base.bed.gz        — per-base depth (optional, --no-per-base to skip)
      *.regions.bed.gz         — per-window mean depth (--by N)

    Notes:
      - mosdepth is substantially faster than samtools depth or bedtools genomecov
        due to its CIGAR-based implementation in Nim.
      - --quantize: bin depth values to reduce output file size (useful for WGS).
      - --by: window size for regional depth (default 200 bp; increase to 1000 for speed).
      - For T2T-CHM13, chrM coverage is often 100–1000× due to mitochondrial DNA.
        This does not affect the nuclear genome coverage calculation when mosdepth
        reports are used with --chrom filtering.
========================================================================================
*/

process MOSDEPTH {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/coverage/${meta.id}", mode: 'copy'

    container   'sktrinh12/mosdepth:latest'

    input:
    tuple val(meta), path(bam), path(bai)
    path  reference_fai
    val   window_size           // window size in bp for regional coverage

    output:
    tuple val(meta), path("${meta.id}.mosdepth.summary.txt"),        emit: summary
    tuple val(meta), path("${meta.id}.mosdepth.global.dist.txt"),    emit: global_dist
    tuple val(meta), path("${meta.id}.regions.bed.gz"),              emit: regions,    optional: true
    tuple val(meta), path("${meta.id}.regions.bed.gz.csi"),          emit: regions_csi, optional: true
    path  "versions.yml",                                             emit: versions

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def win     = window_size ?: 200

    // --no-per-base: skip per-base output (large file, usually not needed for QC)
    // --quantize: merge adjacent positions with same depth bin (reduces file size)
    // --by: window size for regional depth BED
    // Quantize bins: 0:4 = very low, 4:10 = low, 10:50 = normal, 50: = high coverage
    """
    mosdepth \\
        --threads ${task.cpus} \\
        --by ${win} \\
        --no-per-base \\
        --quantize 0:4:10:50: \\
        --fasta ${reference_fai} \\
        ${prefix} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.mosdepth.summary.txt
    touch ${meta.id}.mosdepth.global.dist.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: "stub"
    END_VERSIONS
    """
}
