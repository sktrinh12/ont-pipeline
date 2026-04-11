/*
============================================================================================
    MODULE: ODGI_VIZ
    Tool   : odgi (>=0.9.0)
    Purpose: Visualize pangenome graph structure from ODGI format graph.
             Shows path coverage across the pangenome.
    -s: for colour coding separated string by '#'
    -region_arg: exactly which part of your pangenome graph want to snap a picture of; path_name:start-end
============================================================================================
*/

process ODGI_VIZ {
    tag         "$meta.id"
    label       'process_low'
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
        -s '#' \\
        -P \\
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

/*
============================================================================================
    MODULE: ODGI_BUILD
    Tool   : vg (>=1.73), odgi (>=0.9.0)
    This grep finds any Path (P) line that ends in a '*' (empty node list) 
    and deletes it so ODGI doesn't crash.
    Build the odgi graph
============================================================================================
*/

process ODGI_BUILD {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/pangenome", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(augmented_vg)

    output:
    tuple val(meta), path("${meta.id}.og"), emit: odgi
    path "versions.yml",                  emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vg convert ${augmented_vg} -f -W > ${prefix}.gfa
    grep -vP "P\\t.*\\t\\*" ${prefix}.gfa > ${prefix}.clean.gfa
    odgi build -g ${prefix}.gfa -o ${prefix}.og -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | head -1 | awk '{print \$2}')
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.og
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: "stub"
        odgi: "stub"
    END_VERSIONS
    """
}
