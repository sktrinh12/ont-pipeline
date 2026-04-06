/*
========================================================================================
    MODULE: MULTIQC
    Tool   : MultiQC (>=1.21)
    Purpose: Aggregate QC results from multiple tools into a single HTML report.
========================================================================================
*/

process MULTIQC {
    label       'process_low'
    publishDir  "${params.outdir}/multiqc", mode: 'copy'
    container   'multiqc/multiqc:latest'

    input:
    path qc_files
    val  title

    output:
    path "multiqc_report.html",  emit: report
    path "multiqc_data/",        emit: data,     optional: true
    path "versions.yml",         emit: versions

    script:
    def title_arg = title ? "--title \"${title}\"" : ''
    """
    echo "=== Debug: Files received by MultiQC ==="
    ls -la

    echo ""
    echo "=== Running MultiQC ==="
    multiqc \\
        --force \\
        --ignore-symlinks \\
        ${title_arg} \\
        -n multiqc_report.html \\
        -o . \\
        . 2>&1 | tee multiqc_run.log

    if [ ! -f multiqc_report.html ]; then
        echo "ERROR: MultiQC found no data to process" >&2
        echo "Files in work dir:" >&2
        ls -la >&2
        echo "MultiQC log:" >&2
        cat multiqc_run.log >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_data
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: "stub"
    END_VERSIONS
    """
}
