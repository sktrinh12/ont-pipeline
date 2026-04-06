/*
========================================================================================
    MODULE: BCFTOOLS_NORM
    Purpose: Left-align indels, split multiallelic sites, normalise representation.
             This is required before annotation (VEP, ANNOVAR) and merging VCFs.

    Notes:
      - -m -any: split multiallelic records into biallelic (one alt per line).
        Required for Clair3 output which may produce multi-allelic records.
      - -f: reference FASTA for left-alignment.
      - --check-ref: warn when REF allele doesn't match reference (useful for
        catching reference mismatch errors).
      - Do NOT use bcftools norm on Sniffles2 SV VCFs unless you understand
        the implications — SV representations can be altered by normalisation.
========================================================================================
*/

process BCFTOOLS_NORM {
    tag         "${meta.id} [${mode}]"
    label       'process_low'

    container   'sktrinh12/bcftools:latest'

    input:
    tuple val(meta), path(vcf)
    path  reference
    val   mode          // 'snp_indel' or 'sv' — controls normalisation flags

    output:
    tuple val(meta), path("${meta.id}.${mode}.norm.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.${mode}.norm.vcf.gz.tbi"), emit: tbi
    path  "versions.yml",                                         emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def norm_opt = mode == 'snp_indel'
        ? "-m -any -f ${reference} --check-ref w"
        : "-m -any"   // SVs: split multiallelics but do NOT left-align

    """
    bcftools norm \\
        ${norm_opt} \\
        --threads ${task.cpus} \\
        -O z \\
        -o ${prefix}.${mode}.norm.vcf.gz \\
        ${vcf}

    bcftools index --tbi ${prefix}.${mode}.norm.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${mode}.norm.vcf.gz
    touch ${meta.id}.${mode}.norm.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: BCFTOOLS_FILTER
    Purpose: Apply hard filters to VCF. Retains only PASS variants above QUAL threshold.

    Notes:
      - For Clair3 output, FILTER=PASS is applied by Clair3 internally.
        bcftools filter here removes any records where FILTER != PASS (e.g., LowQual).
      - QUAL >= 20 corresponds to 99% genotype confidence (Phred-scaled).
      - SNPs and Indels may benefit from different QUAL thresholds; Clair3's FILTER
        tags are generally more reliable than raw QUAL for model-based callers.
========================================================================================
*/

process BCFTOOLS_FILTER {
    tag         "${meta.id} [${mode}]"
    label       'process_low'
    publishDir  "${params.outdir}/variants/${mode}/${meta.id}", mode: 'copy'

    container   'sktrinh12/bcftools:latest'

    input:
    tuple val(meta), path(vcf)
    val   filter_expr   // e.g., 'FILTER=="PASS" && QUAL>=20'
    val   mode          // label for output prefix

    output:
    tuple val(meta), path("${meta.id}.${mode}.filtered.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.${mode}.filtered.vcf.gz.tbi"), emit: tbi
    path  "versions.yml",                                             emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view \\
        --include ${filter_expr} \\
        --threads ${task.cpus} \\
        -O z \\
        -o ${prefix}.${mode}.filtered.vcf.gz \\
        ${vcf}

    bcftools index --tbi ${prefix}.${mode}.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${mode}.filtered.vcf.gz
    touch ${meta.id}.${mode}.filtered.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: BCFTOOLS_STATS
    Purpose: Generate variant statistics for MultiQC integration.
             Produces counts of SNPs, Indels, transitions/transversions, etc.
========================================================================================
*/

process BCFTOOLS_STATS {
    tag         "${meta.id} [${mode}]"
    label       'process_low'
    publishDir  "${params.outdir}/qc/variant_stats/${meta.id}", mode: 'copy'

    container   'sktrinh12/bcftools:latest'

    input:
    tuple val(meta), path(vcf)
    val   mode

    output:
    tuple val(meta), path("${meta.id}.${mode}.bcftools_stats.txt"), emit: stats
    path  "versions.yml",                                           emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools stats \\
        --threads ${task.cpus} \\
        ${vcf} \\
        > ${prefix}.${mode}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${mode}.bcftools_stats.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: "stub"
    END_VERSIONS
    """
}
