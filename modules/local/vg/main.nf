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
    publishDir  "${params.outdir}/pangenome/${meta.id}", mode: 'copy',
                saveAs: { fn -> fn.endsWith('.stats') ? fn : null }

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
    """
    vg giraffe \
        -t ${task.cpus} \
        -Z ${gbz} \
        -m ${min_index} \
        -d ${dist_index} \
        -f ${reads} \
        --output-format gam \
        --zipcode-name ${zipcode} \
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
