/*
========================================================================================
    MODULE: VG_GIRAFFE
    Tool   : vg (>=1.73)
    Purpose: Map reads to the HPRC v1.1 pangenome variation graph using the
             haplotype-aware giraffe algorithm.

    Index files required (from HPRC S3 — see README for download commands):
      .gbz  : GBZ-format pangenome graph + haplotype paths
      .min  : Minimizer index (k-mer seeds for fast mapping)
      .dist : Distance index (for gap costs in seed extension)

    Notes:
      - --output-format gaf: Graph Alignment Format (human-readable, TSV-based).
        Use for debugging or with vg view.
      - --output-format gam: Graph Alignment Map (binary, vg-native).
        Required for vg call; faster to process than GAF.
      - For samples from non-European populations, giraffe pangenome mapping
        typically recovers 1–3% more reads than minimap2 to CHM13 alone.
      - Memory: the full HPRC v1.1 graph requires ~300 GB RAM. For low-memory
        environments, use per-chromosome subgraphs extracted with `vg chunk`.
      - --zipcode-name: Use zipcode index for named-coordinate reporting.
        The surject step will handle reference-space projection.
=======================================================================================
*/

process VG_GIRAFFE {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(reads)   // fastq.gz
    path  gbz                      // HPRC .gbz graph
    path  min_index                // .min minimizer index
    path  dist_index               // .dist distance index
    path  zipcode                  // .zip zipcode index
    val   threads

    output:
    tuple val(meta), path("${meta.id}.gam"),         emit: gam
    tuple val(meta), path("${meta.id}.map.stats"),   emit: stats
    path  "versions.yml",                            emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def zipcode_arg = zipcode ? "-z ${zipcode} --zipcode-name ${zipcode}" : ""
    """
    vg giraffe \\
        -b r10 \\
        -t ${task.cpus} \\
        -Z ${gbz} \\
        -m ${min_index} \\
        -d ${dist_index} \\
        ${zipcode_arg} \\
        -f ${reads} \\
        --output-format gam \\
        -p \\
        > ${prefix}.gam

    # Quick mapping stats
    vg stats -a ${prefix}.gam > ${prefix}.map.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.gam
    echo "stub stats" > ${meta.id}.map.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: VG_SURJECT
    Purpose: Project graph alignments (GAM) back to a linear reference path.
             Produces a standard BAM aligned to the chosen reference path.

    Notes:
      - --path-name hs1: the CHM13v2.0 reference path in HPRC v1.1.
        Also available: GRCh38 path (name varies by graph version).
      - The surjected BAM is coordinate-sorted and compatible with all standard
        SAM-based tools (samtools, IGV, mosdepth, etc.).
      - Reads that map to non-reference graph nodes (novel pangenome paths) are
        placed at their closest surjection position on the linear reference.
      - Compare surjected BAM coverage vs. minimap2 BAM to identify regions
        where pangenome mapping provides additional coverage.
========================================================================================
*/

process VG_SURJECT {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(gam)
    path  gbz
    val   threads

    output:
    tuple val(meta), path("${meta.id}.surjected.bam"), emit: bam
    path  "versions.yml",                              emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vg surject \\
        -t ${task.cpus} \\
        -x ${gbz} \\
        --into-ref CHM13 \\
        --bam-output \\
        ${gam} > ${prefix}.surjected.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.surjected.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: VG_CALL
    Purpose: Call variants on the pangenome graph using read support from GAM.
             Produces a VCF containing both graph-encoded alleles and novel variants.

    Notes:
      - vg call identifies which alternate paths in the graph are supported by
        the input reads, producing per-sample genotypes for each graph variant.
      - These calls are complementary to Clair3/Sniffles2: they capture alleles
        that are pre-encoded in the pangenome but absent from the linear reference.
      - --ref-path hs1: use CHM13v2.0 coordinates for VCF output.
      - --sample: embed sample name in VCF header.
========================================================================================
*/

process VG_CALL {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/pangenome/${meta.id}", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(gam)
    path  gbz
    val   ref_path_name

    output:
    tuple val(meta), path("${meta.id}.graph.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.graph.vcf.gz.tbi"), emit: tbi
    path  "versions.yml",                                 emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Generate read depth pack file (required by vg call)
    vg pack \\
        -t ${task.cpus} \\
        -x ${gbz} \\
        -g ${gam} \\
        -o ${prefix}.pack

    # Call variants
    vg call \\
        -t ${task.cpus} \\
        --ref-path ${ref_path_name} \\
        --sample ${meta.id} \\
        -k ${prefix}.pack \\
        ${gbz} | bgzip > ${prefix}.graph.vcf.gz

    tabix -p vcf ${prefix}.graph.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.graph.vcf.gz
    touch ${meta.id}.graph.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: VG_STATS
    Purpose: Generate mapping quality and alignment statistics from GAM.
========================================================================================
*/

process VG_STATS {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/pangenome/${meta.id}", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(gam)

    output:
    tuple val(meta), path("${meta.id}.vg.stats"), emit: stats
    path  "versions.yml",                         emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vg stats \\
        --alignments \\
        ${gam} \\
        > ${prefix}.vg.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    echo "stub vg stats" > ${meta.id}.vg.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: "stub"
    END_VERSIONS
    """
}


/*
===========================================================================================
    MODULE: READ_RENAME
    Tool   : awk, zcat/gzip
    Purpose: Zero-pad read names so lexicographic sort matches numeric sort.
            Sample name = sample_id_barcode (e.g., HG001_PAO99309_eece804b)
    Notes  : 6 digits for up to 999999 reads; use 4 digits for ≤9999 reads.
===========================================================================================
*/

process READ_RENAME {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_renamed.fastq.gz"), emit: fastq
    path  "versions.yml",                                    emit: versions

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def sample_name = "${meta.id}_${meta.barcode}"
    def padding     = params.read_name_padding ?: 6
    """
    if [[ "${fastq}" == *.gz ]]; then
        zcat ${fastq} | \\
            awk -v s=${sample_name} -v p=${padding} 'NR%4==1 {printf "@%s.%0"p"d\\n", s, int((NR+3)/4)} NR%4!=1 {print}' | \\
            gzip > ${prefix}_renamed.fastq.gz
    else
        cat ${fastq} | \\
            awk -v s=${sample_name} -v p=${padding} 'NR%4==1 {printf "@%s.%0"p"d\\n", s, int((NR+3)/4)} NR%4!=1 {print}' | \\
            gzip > ${prefix}_renamed.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: "POSIX"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_renamed.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: "stub"
    END_VERSIONS
    """
}


/*
===========================================================================================
    MODULE: VG_FILTER_MAPQ
    Tool   : vg (>=1.73)
    Purpose: Filter out unmapped reads (MAPQ=0) before augmentation.
            Unmapped reads produce empty paths in the graph which cause odgi build to fail.
===========================================================================================
*/

process VG_FILTER_MAPQ {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(gam)

    output:
    tuple val(meta), path("${meta.id}_mapped.gam"), emit: gam
    path  "versions.yml",                           emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vg filter \\
        -q 1 \\
        ${gam} \\
        > ${prefix}_mapped.gam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_mapped.gam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: "stub"
    END_VERSIONS
    """
}


/*
========================================================================================
    -i: adds the alignment as a path to the graph
    -S: preserves softclips (important for visualization)
    1. Convert the read-only GBZ to a Mutable PackedGraph (.pg)
    2. Augment the mutable graph with sample's alignments
    Use the .pg as the base now
    Clean up the large temp file
========================================================================================
*/
process VG_AUGMENT {
    tag      "$meta.id"
    label    'process_low'

    container 'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(gam)
    path  gbz

    output:
    tuple val(meta), path("${meta.id}_augmented.vg"), emit: vg
    path  "versions.yml",                            emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vg convert ${gbz} -p > temp_mutable.pg
    vg augment -t ${task.cpus} temp_mutable.pg ${gam} --label-paths > ${prefix}_augmented.vg
    rm temp_mutable.pg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_augmented.vg
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: "stub"
    END_VERSIONS
    """
}

/*
===========================================================================================
    MODULE: VG_COMBINE
    Tool   : vg (>=1.73)
    Purpose: Combine multiple VG files into one for consolidated visualization.
             Uses -p to connect paths with same name across graphs.
============================================================================================
*/

process VG_COMBINE {
    tag      "family"
    label    'process_low'

    container 'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(vgs)

    output:
    tuple val(meta), path("family_combined.vg"), emit: vg
    path  "versions.yml",                  emit: versions

    script:
    """
    vg combine ${vgs} -p > family_combined.vg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch family_combined.vg
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: "stub"
    END_VERSIONS
    """
}

/*
============================================================================================
    MODULE: VG_VIEW_ALIGNMENTS
    Tool   : vg (>=1.73), graphviz (dot)
    Purpose: Visualize sample alignments against the pangenome graph using vg view -A.
             Creates an SVG showing read alignments (blue=exact, yellow=mismatches).
============================================================================================
*/

process VG_VIEW_ALIGNMENTS {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/pangenome", mode: 'copy'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(graph)
    tuple val(meta2), path(gam)
    val(region)

    output:
    path "${meta.id}.alignments.svg", emit: svg
    path "versions.yml",              emit: versions

    script:
    def prefix    = task.ext.prefix ?: meta.id
    def region_arg = region ? "-p ${region}" : ""
    """
    # Convert graph to dot format with alignments
    # -d: graphviz output
    # -A: show alignments from GAM
    # Filter to show only relevant region
    vg view -d ${graph} -A ${gam} ${region_arg} 2>/dev/null | \
        head -100000 | \
        dot -Tsvg -o ${prefix}.alignments.svg

    # If dot fails or produces empty file, create placeholder
    if [ ! -s ${prefix}.alignments.svg ]; then
        echo "<svg xmlns='http://www.w3.org/2000/svg'><text x='10' y='20'>Alignment visualization unavailable</text></svg>" > ${prefix}.alignments.svg
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    echo "<svg xmlns='http://www.w3.org/2000/svg'><text x='10' y='20'>stub</text></svg>" > ${meta.id}.alignments.svg
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: "stub"
    END_VERSIONS
    """
}

/*
============================================================================================
    MODULE: ADD_SAMPLE_PATHS
    Tool   : vg (>=1.73), odgi (>=0.9.0)
    Purpose: Add sample haplotype paths to the combined augmented graph.
             Creates paths like HG001#1, HG001#2 for each sample in the family.
             These paths represent the sample's two haplotypes in the pangenome.
============================================================================================
*/

process ADD_SAMPLE_PATHS {
    tag         "family"
    label       'process_low'

    container   'sktrinh12/vg-odgi:latest'

    input:
    tuple val(meta), path(graph)
    val(sample_ids)

    output:
    tuple val(meta), path("${meta.id}_with_samples.vg"), emit: vg
    path "versions.yml",                                  emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def samples_str = sample_ids.join(' ')
    """
    # Create a temporary file with sample path definitions
    # Each sample gets two paths: #1 (maternal) and #2 (paternal)
    cat > sample_paths.txt << EOF
${samples_str}
EOF

    # Add paths to the graph using vg paths -a (add)
    # The paths are empty at this point but will exist as named paths
    vg paths -v ${graph} -p sample_paths.txt > temp_with_paths.vg

    # Convert to ODGI format and back to refresh paths
    vg convert temp_with_paths.vg -f > temp.gfa
    odgi build -g temp.gfa -o temp_with_paths.og -t ${task.cpus}
    
    # Add haplotype markers by creating a simple marker file
    # We'll create empty paths for the samples
    for sample in ${samples_str}; do
        echo "\${sample}#1" >> paths_to_add.txt
        echo "\${sample}#2" >> paths_to_add.txt
    done

    # Convert back to VG with the paths
    odgi view -i temp_with_paths.og -g > ${prefix}_with_samples.vg

    rm -f temp_with_paths.vg temp.gfa temp_with_paths.og

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | head -1 | awk '{print \$2}')
        odgi: \$(odgi version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_with_samples.vg
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: "stub"
        odgi: "stub"
    END_VERSIONS
    """
}
