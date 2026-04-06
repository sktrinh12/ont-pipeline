/*
========================================================================================
    MODULE: CLAIR3
    Tool   : Clair3 (>=1.0.4)
    Purpose: Call SNPs and Indels from ONT long-read alignments using deep learning.
    Notes  :
      - Clair3 requires a model trained for the specific chemistry + accuracy tier.
        Models are bundled in the Docker image or can be downloaded from:
        https://github.com/HKU-BAL/Clair3#pre-trained-models
      - --gvcf: enables GVCF output for cohort joint genotyping via GLnexus.
      - --platform ont: forces ONT-specific calling parameters.
      - Clair3 runs chromosome-parallel internally; the --threads flag controls
        per-chromosome parallelism within a single process.
========================================================================================
*/

process CLAIR3 {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/variants/snp_indel/${meta.id}", mode: 'copy'

    container   'hkubal/clair3:latest'

    input:
    tuple val(meta), path(bam), path(bai)
    path  reference
    path  reference_fai
    val   model_path
    val   min_coverage

    output:
    tuple val(meta), path("clair3_${meta.id}/merge_output.vcf.gz"),     emit: vcf
    tuple val(meta), path("clair3_${meta.id}/merge_output.gvcf.gz"),    emit: gvcf,    optional: true
    tuple val(meta), path("clair3_${meta.id}/"),                        emit: output_dir
    path  "versions.yml",                                                emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def gvcf_opt = '--gvcf'   // always produce GVCF for cohort compatibility

    // --model_path: model must match basecalling chemistry.
    //   Clair3 ships models in /opt/models/ inside the container.
    //   Bundled models: ont_sup_g5, ont_hac_g5, ont_fast_g5
    //   For newer R10.4.1 v4.3 data, download the appropriate model from the GitHub.
    //
    // --include_all_ctgs: call on all contigs, not just chr1-22,X,Y.
    //   Important for T2T-CHM13 which includes chrM and unplaced scaffolds.
    """
    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --platform=ont \\
        --model_path=/opt/models/${model_path} \\
        --output=clair3_${prefix} \\
        --min_coverage=${min_coverage} \\
        --include_all_ctgs \\
        --no_phasing_for_fa \\
        ${gvcf_opt}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p clair3_${meta.id}
    touch clair3_${meta.id}/merge_output.vcf.gz
    touch clair3_${meta.id}/merge_output.gvcf.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: "stub"
    END_VERSIONS
    """
}
