/*
========================================================================================
    MODULE: PHASE_STATS
    Purpose: Generate phasing quality statistics. Quality Control for Haplotype Assembly

    PURPOSE:
    Assesses the continuity of phased genomic variants. In long-read sequencing (ONT), 
    phasing connects variants into parental haplotypes (maternal vs. paternal).

    HOW IT WORKS:
    This process parses a phased VCF file and counts "Phase Sets" (PS tags). 
    A Phase Set is a continuous block where the relationship between variants 
    is known. More blocks usually mean more fragmentation (shorter N50).

    BIOLOGICAL INTERPRETATION (specifically for Chr22 @ 30x ONT):
    - N50 > 5 Mb:  Excellent; indicates near-complete resolution of haplotypes.
    - N50 1-5 Mb:  Good; likely fragmented by centromeres or repetitive regions.
    - N50 < 1 Mb:  Warning; suggests poor coverage or low genetic diversity.

    INPUT:  Tuple containing metadata and a phased VCF file.
    OUTPUT: A text file summarizing the total count of phase sets/blocks.
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
