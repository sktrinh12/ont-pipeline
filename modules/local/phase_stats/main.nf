/*
========================================================================================
    MODULE: PHASE_STATS
    Purpose: Generate phasing quality statistics.
========================================================================================
*/

process PHASE_STATS {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/phasing/${meta.id}", mode: 'copy'

    container   'sktrinh12/bcftools:latest'

    input:
    tuple val(meta), path(phased_vcf)

    output:
    tuple val(meta), path("${meta.id}.phase_stats.txt"), emit: stats
    path  "versions.yml",                                emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view -H ${phased_vcf} | awk -F'\\t' 'BEGIN{ps=0; blocks=0} {if(\$8 ~ /PS=/) {ps++; if(\$8 ~ /PS=([0-9]+)/) blocks++}} END{print "Phase sets: " ps "\\nPhase blocks: " blocks}' > ${prefix}.phase_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.phase_stats.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
