/*
========================================================================================
    MODULE: LONGPHASE_PHASE
    Tool   : longphase (>=1.6)
    Purpose: Phase heterozygous SNPs, Indels, and optionally SVs from an aligned
             ONT BAM and a Clair3 VCF. Produces a VCF with PS (phase set) tags
             encoding haplotype block assignments.

    Algorithm summary:
      longphase uses a read-backed phasing approach:
        1. For each heterozygous variant, reads spanning multiple het sites are
           used to connect the variants into phase blocks (PS).
        2. With long ONT reads (N50 >10 kb), individual reads routinely span
           dozens of het sites, producing very long phase blocks (N50 >10 Mb).
        3. SVs can be included via --sv-file; junction-spanning reads link SVs
           to SNPs, dramatically improving phase block continuity in SV-dense
           genomic regions (e.g., segmental duplications).

    Notes:
      - Input VCF must be indexed (.tbi or .csi).
      - The --ont flag activates ONT-specific parameters (different from PacBio HiFi).
      - --indels: include indel phasing (default is SNPs only; indels add noise
        but also improve block connectivity for high-coverage samples).
      - Output: phased VCF with PS (phase set), GT with pipe notation (0|1 vs 0/1).
========================================================================================
*/

process LONGPHASE_PHASE {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/longphase:latest'

    input:
    tuple val(meta), path(bam), path(bai), path(vcf)
    path  reference
    path  reference_fai
    val   threads

    output:
    tuple val(meta), path("${meta.id}.phased.vcf.gz"),     emit: phased_vcf
    tuple val(meta), path("${meta.id}.phased.vcf.gz.tbi"), emit: tbi
    path  "versions.yml",                                  emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    longphase phase \\
        --ont \\
        --indels \\
        -s ${vcf} \\
        -b ${bam} \\
        -r ${reference} \\
        -t ${threads} \\
        -o ${prefix}.phased

    # Compress and index output
    /opt/conda/bin/bgzip -@ ${task.cpus} ${prefix}.phased.vcf
    /opt/conda/bin/tabix -p vcf ${prefix}.phased.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version 2>&1 | awk '{print \$2}')
        bgzip:     \$(/opt/conda/bin/bgzip --version 2>&1 | head -1 | awk '{print \$2}')
        tabix:     \$(/opt/conda/bin/tabix --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.phased.vcf.gz
    touch ${meta.id}.phased.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: LONGPHASE_HAPLOTAG
    Purpose: Physical grouping of sequencing reads into parental haplotypes.
             Assigns 'HP' tags (HP:1, HP:2) to reads by matching them to a phased VCF.

    BIOLOGICAL UTILITY:
      - Allele-Specific Methylation (ASM): Enables 'modkit --partition-tag HP' to 
        detect differential methylation at imprinted or X-inactivated loci.
      - Phased Variant Calling: Improves structural variant (SV) detection.
      - IGV Visualization: Allows "Group by: HP" to see distinct parental read streams.

    EXPECTED PERFORMANCE (ONT 30x R10.4.1):
      - ~85–92% Tagging Rate: Standard for human WGS.
      - HP:0 (Unassigned): Occurs in homozygous regions or low-complexity repeats 
        where a lack of heterozygous SNPs prevents parental assignment.
========================================================================================
*/

process LONGPHASE_HAPLOTAG {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/phasing/${meta.id}", mode: 'copy'

    container   'sktrinh12/longphase:latest'

    input:
    tuple val(meta), path(bam), path(bai), path(phased_vcf)
    path  reference
    path  reference_fai

    output:
    tuple val(meta), path("${meta.id}.haplotagged.bam"), emit: bam
    path  "versions.yml",                                emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    longphase haplotag \\
        -s ${phased_vcf} \\
        -b ${bam} \\
        -r ${reference} \\
        -t ${task.cpus} \\
        -o ${prefix}.haplotagged

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.haplotagged.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: "stub"
    END_VERSIONS
    """
}
