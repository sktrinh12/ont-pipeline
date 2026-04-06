/*
============================================================================================
    MODULE: ODGI_VIZ
    Tool   : odgi (>=0.9.0)
    Purpose: Visualize pangenome graph structure from ODGI format graph.
             Shows path coverage across the pangenome.
============================================================================================
*/

process ODGI_VIZ {
    tag         "$meta.id"
    label       'process_medium'
    publishDir  "${params.outdir}/qc/pangenome", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(odgi)
    val(region)

    output:
    path "${meta.id}.viz.png", emit: png
    path "versions.yml",      emit: versions

    script:
    def prefix     = task.ext.prefix ?: meta.id
    def region_arg = region ? "-r ${region}" : ""
    """
    odgi viz \\
        -i ${odgi} \\
        -o ${prefix}.viz.png \\
        ${region_arg} \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.viz.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}
