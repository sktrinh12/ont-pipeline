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
    Purpose: Convert GBZ pangenome graph to ODGI format for visualization.
    Direct conversion from GBZ to ODGI
============================================================================================
 TODO: Enable alignment overlay
       1. Add GAM to input tuple: tuple val(meta), path(gbz), path(gam)
       2. Add vg augment + vg convert steps below
       3. Update pangenome.nf to pass GAM to ODGI_BUILD

       - Issue of Computational Scale.
       - If "augment" a graph with several gigabytes of FASTQ data, .odgi file will become massive. Instead of having ~50 paths (the HPRC references), it will have millions of paths (one for every read).

       - Don't use odgi viz -A for the whole trio across all of Chr22. It will likely crash nodes or produce a PNG that is 100,000 pixels tall.

       - First, run current odgi_region = 'CHM13#chr22:0-*' to see the base graph.
       - Then, if find a specific region (like a 10kb window where you suspect a structural variant), use a smaller odgi_region and then try the vg augment path.
       - Need to add: Augment (The missing link)
         // This turns the .gam alignments into new paths in the graph
         VG_AUGMENT(VG_GIRAFFE.out.gam, ch_gbz)
============================================================================================
*/

process ODGI_BUILD {
    tag         "$meta.id"
    label       'process_low'
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
    vg convert -g ${gbz} -o > ${prefix}.og

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
