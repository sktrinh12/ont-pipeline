/*
=======================================================================================
    MODULE: SV_SUMMARY
    Purpose: Generate summary statistics for structural variants.
=======================================================================================
*/

process SV_SUMMARY {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/sv_stats/${meta.id}", mode: 'copy'

    container   'sktrinh12/python:3.12'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.sv_summary.txt"), emit: summary
    path  "versions.yml",                               emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_gzipped = vcf.toString().endsWith('.gz')
    """
    python3 -c '
import gzip
vcf = "${vcf}"
counts = {"DEL": 0, "INS": 0, "INV": 0, "DUP": 0, "BND": 0}
${is_gzipped ? "with gzip.open(vcf, \"rt\") as f:" : "with open(vcf) as f:"}
    for line in f:
        if line.startswith("#"): continue
        fields = line.split("\t")
        info = fields[7]
        for svtype in counts:
            if "SVTYPE=" + svtype in info:
                counts[svtype] += 1
                break
with open("${prefix}.sv_summary.txt", "w") as out:
    out.write("SV Type Counts:\\n")
    for t, c in counts.items():
        out.write(f"  {t}: {c}\\n")
'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sv_summary.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
