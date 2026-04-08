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
        -s '#' \\
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
    Purpose: Convert GBZ pangenome graph to ODGI format for visualization.
             1. GBZ → GFA using vg view
             2. GFA → ODGI format using odgi build
============================================================================================
 TODO: Enable alignment overlay (requires VG issue #4240 fix)
       To re-enable:
       1. Add GAM to input tuple: tuple val(meta), path(gbz), path(gam)
       2. Uncomment vg augment + vg convert steps below
       3. Update pangenome.nf to pass GAM to ODGI_BUILD
============================================================================================
*/

process ODGI_BUILD {
    tag         "$meta.id"
    label       'process_medium'
    publishDir  "${params.outdir}/qc/pangenome", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(gbz)

    output:
    tuple val(meta), path("${meta.id}.og"), emit: odgi
    path "versions.yml",                  emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vg view --gfa ${gbz} > ${prefix}.gfa

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
