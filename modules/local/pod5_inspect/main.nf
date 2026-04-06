/*
========================================================================================
    MODULE: POD5_INSPECT
    Tool   : pod5 (>=0.1.0)
    Purpose: Inspect POD5 files, extract read count, signal stats, detect truncated files
========================================================================================
*/

process POD5_INSPECT {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/pod5_inspect/${meta.id}", mode: 'copy'

    container   'community.wave.seqera.io/library/pod5:latest'

    input:
    tuple val(meta), path(pod5_dir)

    output:
    tuple val(meta), path("${meta.id}_pod5_summary.txt"), emit: summary
    path  "versions.yml",                              emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pod5 inspect ${pod5_dir} > ${prefix}_pod5_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pod5: \$(pod5 --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_pod5_summary.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pod5: "stub"
    END_VERSIONS
    """
}
