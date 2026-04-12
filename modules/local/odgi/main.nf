/*
===========================================================================================
    MODULE: ODGI_VIZ_HEATMAP
    Tool   : odgi (>=0.9.0)
    Purpose: Visualize pangenome graph - Compressed Coverage Heatmap.
             Collapses all read paths into a single heatmap row alongside reference haplotypes.
             Shows where reads covered the pangenome.
    Flags  : -O = compressed mode (cannot combine with -p)
             -I = ignore prefix to exclude (e.g., GRCh38)
             -P = pixel-based rendering
===========================================================================================
*/

process ODGI_VIZ_HEATMAP {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/pangenome", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(odgi)
    val(region)
    val(ignore_prefix)

    output:
    path "${meta.id}.heatmap.png", emit: heatmap_png
    path "versions.yml",           emit: versions

    script:
    def prefix       = task.ext.prefix ?: meta.id
    def region_arg  = region && region != 'null' ? "-r ${region}" : ""
    def ignore_arg  = ignore_prefix && ignore_prefix != 'null' ? "-I ${ignore_prefix}" : ""
    """
    odgi viz \
        -i ${odgi} \
        -o ${prefix}.heatmap.png \
        ${region_arg} \
        ${ignore_arg} \
        -O \
        -P \
        -x 2000 -y 300 \
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.heatmap.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}


/*
===========================================================================================
    MODULE: ODGI_VIZ_DEPTH
    Tool   : odgi (>=0.9.0)
    Purpose: Visualize pangenome graph - Per-Read Depth Coloring.
             Colors bins black→blue by mean depth. Shows coverage uniformity.
             Uses display path list to restrict to CHM13 + sample reads only.
    Flags  : -m = mean depth coloring
             -p = path list file
             -P = pixel-based rendering
===========================================================================================
*/

process ODGI_VIZ_DEPTH {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/pangenome", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(odgi)
    path(path_list)
    val(region)

    output:
    path "${meta.id}.depth.png", emit: depth_png
    path "versions.yml",         emit: versions

    script:
    def prefix      = task.ext.prefix ?: meta.id
    def region_arg  = region && region != 'null' ? "-r ${region}" : ""
    """
    odgi viz \
        -i ${odgi} \
        -o ${prefix}.depth.png \
        -p ${path_list} \
        ${region_arg} \
        -m \
        -P \
        -x 2000 -y 500 \
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.depth.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}


/*
===========================================================================================
    MODULE: ODGI_VIZ_CLUSTERED
    Tool   : odgi (>=0.9.0)
    Purpose: Visualize pangenome graph - Clustered View.
             Orders reads by pairwise Jaccard similarity so reads covering the same regions
             cluster together. Most structured view — reveals coverage patterns across haplotypes.
    Flags  : -k = cluster by Jaccard similarity (cannot combine with -p)
             -I = ignore prefix to exclude (e.g., GRCh38)
             -P = pixel-based rendering
===========================================================================================
*/

process ODGI_VIZ_CLUSTERED {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/pangenome", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(odgi)
    val(region)
    val(ignore_prefix)

    output:
    path "${meta.id}.clustered.png", emit: clustered_png
    path "versions.yml",              emit: versions

    script:
    def prefix       = task.ext.prefix ?: meta.id
    def region_arg  = region && region != 'null' ? "-r ${region}" : ""
    def ignore_arg  = ignore_prefix && ignore_prefix != 'null' ? "-I ${ignore_prefix}" : ""
    """
    odgi viz \
        -i ${odgi} \
        -o ${prefix}.clustered.png \
        ${region_arg} \
        ${ignore_arg} \
        -k \
        -P \
        -x 2000 -y 1000 \
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.clustered.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}


/*
===========================================================================================
    MODULE: ODGI_BUILD
    Tool   : vg (>=1.73), odgi (>=0.9.0)
    Purpose: Build odgi graph from augmented VG.
    Notes  : - grep removes empty paths (P lines ending in '*') that crash odgi
            - -W forces P-lines (not W-lines) for odgi compatibility
            - -P preserves path information
===========================================================================================
*/

process ODGI_BUILD {
    tag         "$meta.id"
    label       'process_low'

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

    # Safety check: confirm no empty paths (would cause odgi build to fail)
    empty_paths=\$(awk '\$1=="P" && \$3==""' ${prefix}.gfa | wc -l)
    if [ \${empty_paths} -gt 0 ]; then
        echo "ERROR: Found \${empty_paths} empty paths in GFA" >&2
        awk '\$1=="P" && \$3==""' ${prefix}.gfa | head -5
        exit 1
    fi

    # Build odgi graph
    odgi build -g ${prefix}.gfa -o ${prefix}.og -t ${task.cpus} -P

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


/*
===========================================================================================
    MODULE: ODGI_PATHS_FILTER
    Tool   : odgi (>=0.9.0), grep, awk, tr
    Purpose: Create display path list restricting to CHM13 reference + sample reads.
             Used for depth and strand visualizations.
===========================================================================================
*/

process ODGI_PATHS_FILTER {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(odgi)
    val(sample_names)

    output:
    path "display_paths.txt", emit: path_list
    path "versions.yml",      emit: versions

    script:
    def sample_str = sample_names.join('|')
    """
    odgi paths -i ${odgi} -L | \\
        grep -E "^CHM13|${sample_str}" | \\
        tr -d '\\r' | \\
        awk '!seen[\$0]++' > display_paths.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch display_paths.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}