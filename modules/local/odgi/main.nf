/*
==========================================================================================
    MODULE: ODGI_BUILD
    Tool   : vg (>=1.73), odgi (>=0.9.0)
    Purpose: Build odgi graph from augmented VG.
    Notes  : - grep removes empty paths (P lines ending in '*') that crash odgi
            - -W forces P-lines (not W-lines) for odgi compatibility
            - -P preserves path information
==========================================================================================
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
==========================================================================================
    MODULE: ODGI_PATHS_FILTER
    Tool   : odgi (>=0.9.0), grep, awk, tr
    Purpose: Create display path list restricting to CHM13 reference + sample reads.
             Used for depth and strand visualizations.
==========================================================================================
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


/*
==========================================================================================
    MODULE: ODGI_RETAIN
    Tool   : odgi (>=0.9.0)
    Purpose: Extract only specified paths from the full graph to create a smaller,
             focused graph for faster sorting, layout, and visualization.
             Removes HPRC haplotype paths that clutter visualizations.
==========================================================================================
*/

process ODGI_RETAIN {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(full_odgi)
    path(path_list)
    val(region)

    output:
    tuple val(meta), path("${meta.id}.retained.og"), emit: retained_odgi
    path "display_paths_retained.txt",               emit: path_list
    path "versions.yml",                             emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Extract the reference path name from region (e.g. "CHM13#0#chr22" from "CHM13#0#chr22:0-*")
    def ref_path = region ? region.split(':')[0] : ""
    def region_arg = region ? "-r ${region}" : ""
    """
    odgi extract \\
        -i ${full_odgi} \\
        ${region_arg} \\
        -p ${path_list} \\
        -O \\
        -o ${prefix}.retained.og \\
        -P

    # Regenerate display path list from extracted graph
    # (path names gain coordinate suffixes after extraction)
    odgi paths -i ${prefix}.retained.og -L > display_paths_retained.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.retained.og display_paths_retained.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}


/*
==========================================================================================
    MODULE: ODGI_SORT
    Tool   : odgi (>=0.9.0)
    Purpose: Sort graph using PG-SGD algorithm for clean 1D/2D layout.
             Optimizes node ordering to reduce visual complexity.
==========================================================================================
*/

process ODGI_SORT {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(odgi)

    output:
    tuple val(meta), path("${meta.id}.sorted.og"), emit: sorted_odgi
    path "versions.yml",                           emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    odgi sort \\
        -i ${odgi} \\
        -o ${prefix}.sorted.og \\
        -t ${task.cpus} \\
        -P \\
        -p sY \\
        -O
    
    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.og
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}


/*
==========================================================================================
    MODULE: ODGI_LAYOUT
    Tool   : odgi (>=0.9.0)
    Purpose: Compute 2D layout coordinates for graph visualization.
             Uses PG-SGD algorithm to position nodes in 2D space.
==========================================================================================
*/

process ODGI_LAYOUT {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(retained_odgi)

    output:
    tuple val(meta), path("${meta.id}.retained.og.lay"), emit: layout
    path "versions.yml",                               emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    odgi layout \\
        -i ${retained_odgi} \\
        -o ${prefix}.retained.og.lay \\
        -t ${task.cpus} \\
        -P

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.retained.og.lay
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}


/*
==========================================================================================
    MODULE: ODGI_DRAW
    Tool   : odgi (>=0.9.0)
    Purpose: Render 2D graph topology visualization from layout coordinates.
             Shows graph structure: bubbles for variants, tangles for complex regions.
==========================================================================================
*/

process ODGI_DRAW {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/pangenome", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(retained_odgi)
    tuple val(meta2), path(layout)
    val(region)

    output:
    path "${meta.id}.topology.png", emit: topology_png
    path "versions.yml",            emit: versions

    script:
    def prefix      = task.ext.prefix ?: meta.id
    def region_arg  = region && region != 'null' ? "-r ${region}" : ""
    """
    odgi draw \\
        -i ${retained_odgi} \\
        -c ${layout} \\
        -p ${prefix}.topology.png \\
        -C \\
        -H 1000 \\
        -w 1 \\
        -S 2 \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.topology.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}


/*
========================================================================================
    MODULE: ODGI_COMPACT_IDS
    Tool   : odgi (>=0.9.0)
    Purpose: Aggressively compact node IDs to prevent assertion failures in visualization.
             Multiple optimization passes ensure small, contiguous node identifiers.
========================================================================================
*/

process ODGI_COMPACT_IDS {
    tag         '$meta.id'
    label       'process_low'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(odgi)

    output:
    tuple val(meta), path("${meta.id}.compacted.og"), emit: compacted_odgi
    path "versions.yml",                         emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    odgi sort -i ${odgi} -o ${prefix}.compacted.og -p sY -O
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.compacted.og
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}


/*
========================================================================================
    MODULE: ODGI_VIZ_SAMPLE_COVERAGE
    Tool   : odgi (>=0.9.0)
    Purpose: Visualize pangenome graph - Sample-level coverage heatmap.
             Collapses all reads from each sample into a single row using prefix merges.
             Shows which genomic regions are covered by each sample.
    Flags  : -M = prefix merge file (collapses paths with matching prefixes)
              -O = compressed mode (heatmap coloring)
              -P = pixel-based rendering
========================================================================================
*/

process ODGI_VIZ_SAMPLE_COVERAGE {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/pangenome", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(odgi)
    path(path_list)
    path(merge_prefixes)

    output:
    path "${meta.id}.sample_coverage.png", emit: sample_coverage_png
    path "versions.yml",                   emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    odgi viz \\
        -i ${odgi} \\
        -o ${prefix}.sample_coverage.png \\
        -p ${path_list} \\
        -M ${merge_prefixes} \\
        -P \\
        -x 2000 -y 300 \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sample_coverage.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: "stub"
    END_VERSIONS
    """
}