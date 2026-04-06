/*
========================================================================================
    MODULE: CHOPPER_FILTER
    Tool   : Chopper (>=0.7.0)
    Purpose: Quality and length filter FASTQ reads.

    WHY CHOPPER instead of NanoFilt:
      - NanoFilt (https://github.com/wdecoster/nanofilt) is no longer maintained.
        The author (Wouter De Coster) recommends migrating to Chopper.
      - Chopper is a Rust rewrite: 5–20× faster, multi-threaded, lower memory.
      - Chopper reads from stdin/stdout, compatible with pipe-based workflows.
      - NanoFilt had known issues with read header preservation in some edge cases;
        Chopper preserves complete FASTQ headers.

    WHY filter at all?
      - Very short reads (<500 bp) are often chimeric artefacts or adapter dimers.
      - Reads below Q10 (~90% accuracy) increase noise in variant calling,
        particularly for SNP calling in repetitive regions.
      - Filtering BEFORE alignment reduces compute time for all downstream steps.
      - Note: Dorado's --min-qscore already filters at the basecalling stage when
        using pod5 input. Chopper is applied to pre-basecalled FASTQ inputs.
========================================================================================
*/

process CHOPPER_FILTER {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/filtered_reads/${meta.id}", mode: 'copy'

    container   'olistr12/chopper:0.7.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filtered.fastq.gz"), emit: reads
    path  "versions.yml",                                  emit: versions

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def min_len    = params.min_read_length
    def min_qual   = params.min_read_quality
    def max_len    = params.max_read_length > 0 ? "--maxlength ${params.max_read_length}" : ''
    def decompress = reads.toString().endsWith('.gz') ? 'zcat' : 'cat'

    // Pipeline: decompress → chopper filter → gzip output
    // --threads: chopper is multi-threaded (unlike NanoFilt which was single-threaded)
    // --headcrop / --tailcrop: optionally trim adapter bases from read ends
    //   (usually not needed if Dorado trimming is enabled at basecalling)
    """
    ${decompress} ${reads} \\
        | chopper \\
            --quality ${min_qual} \\
            --minlength ${min_len} \\
            ${max_len} \\
            --threads ${task.cpus} \\
        | gzip > ${prefix}.filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.filtered.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "stub"
    END_VERSIONS
    """
}
