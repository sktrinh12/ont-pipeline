/*
========================================================================================
    MODULE: PYCOQC
    Tool   : pycoQC (>=2.5.0)
    Purpose: Generate QC plots from sequencing_summary.txt produced by MinKNOW/Dorado.
             Provides per-channel and per-read-length QC plots.
========================================================================================
*/

process PYCOQC {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/pycoqc/${meta.id}", mode: 'copy'

    container   'community.wave.seqera.io/library/pycoqc:latest'

    input:
    tuple val(meta), path(fast5_dir), path(summary)

    output:
    tuple val(meta), path("${meta.id}_pycoqc.html"), emit: report
    path  "versions.yml",                           emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pycoQC \\
        --summary_file ${summary} \\
        --output_file ${prefix}_pycoqc.html \\
        --min_reads 1000

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pycoqc: \$(pycoQC --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_pycoqc.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pycoqc: "stub"
    END_VERSIONS
    """
}
