/*
========================================================================================
    MODULE: SNIFFLES2
    Tool   : Sniffles2 (>=2.7.3)
    Purpose: Call structural variants (SVs) from ONT long-read alignments.
    SV types supported: DEL, INS, INV, DUP, BND (breakend / translocation)

    Notes  :
      - --tandem-repeats: STRONGLY RECOMMENDED. Without this, Sniffles2 will call
        thousands of false-positive insertions/deletions in tandem repeat regions
        (microsatellites, satellite DNA, centromeres). The BED file must be derived
        from the same reference as the input BAM — use CHM13-specific TRF output.
        Generate with: trf <reference.fa> 2 7 7 80 10 50 500 -f -d -m -ngs
        then parse the .dat output into BED format.

      - --snf output: Sniffles2 Numeric Format — a binary file that encodes all
        read-level SV evidence. Use --snf for single-sample runs, then combine
        .snf files with 'sniffles --input *.snf' for multi-sample joint genotyping.
        This is far superior to merging per-sample VCFs with SURVIVOR because it
        retains genotype-level evidence.

      - --mosaic: adds sensitivity for somatic/mosaic SVs at low allele frequency.
        Use for cancer or developmental mosaicism studies. Increases false-positive
        rate on germline samples.

      - --long-ins-length: default 2500 bp. Long insertions (e.g., TEs, VNTRs) may
        be missed if this is too low. For T2T references with complete centromeres,
        increase to 10000 for maximum sensitivity.

      - CCS/HiFi compatibility: Sniffles2 also works with PacBio HiFi reads.
        Set --minsvlen accordingly (HiFi calls are typically cleaner at ≥50 bp).
========================================================================================
*/

process SNIFFLES2 {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/variants/sv/${meta.id}", mode: 'copy'

    container   'sktrinh12/sniffles2:latest'

    input:
    tuple val(meta), path(bam), path(bai)
    path  reference
    path  tandem_repeats          // BED file (may be empty list [])
    val   min_support
    val   min_svlen

    output:
    tuple val(meta), path("${meta.id}.sv.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.snf"),           emit: snf
    path  "versions.yml",                              emit: versions

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def trf_opt = tandem_repeats ? "--tandem-repeats ${tandem_repeats}" : ''

    // --minsvlen 50: standard minimum SV length; reduces noise from alignment artifacts
    // --minsupport: minimum number of reads supporting SV call
    // --output-rnames: include supporting read names in VCF INFO field (traceability)
    // --long-ins-length 10000: increase sensitivity for large insertions on T2T
    // --reference: required for INS sequence representation in VCF
    """
    sniffles \\
        --input ${bam} \\
        --vcf ${prefix}.sv.vcf.gz \\
        --snf ${prefix}.snf \\
        --reference ${reference} \\
        --threads ${task.cpus} \\
        --minsupport ${min_support} \\
        --minsvlen ${min_svlen} \\
        --long-ins-length 10000 \\
        --output-rnames \\
        ${trf_opt}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: \$(sniffles --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sv.vcf.gz
    touch ${meta.id}.snf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sniffles: "stub"
    END_VERSIONS
    """
}
