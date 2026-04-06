/*
========================================================================================
    MODULE: BEDTOOLS_INTERSECT
    Purpose: Intersect genomic intervals (e.g., bedMethyl with CpG islands).
========================================================================================
*/

process BEDTOOLS_INTERSECT {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/methylation/${meta.id}/annotated", mode: 'copy'

    container   'community.wave.seqera.io/library/bedtools:latest'

    input:
    tuple val(meta), path(a_bed)   // primary BED (e.g., bedMethyl.gz)
    path  b_bed                    // annotation BED (e.g., CpG islands)
    val   flags                    // e.g., '-wa -wb'

    output:
    tuple val(meta), path("${meta.id}.annotated.bed"), emit: result
    path  "versions.yml",                              emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools intersect \\
        -a ${a_bed} \\
        -b ${b_bed} \\
        ${flags} \\
        > ${prefix}.annotated.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.annotated.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: "stub"
    END_VERSIONS
    """
}
