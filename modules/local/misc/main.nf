/*
============================================================================================
    MODULE: MODKIT_SUMMARY
============================================================================================
*/

process MODKIT_SUMMARY {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/methylation/${meta.id}", mode: 'copy'

    container   'ontresearch/modkit:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_modkit_summary.tsv"), emit: summary
    path  "versions.yml",                                   emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    modkit summary \\
        --threads ${task.cpus} \\
        ${bam} \\
        > ${prefix}_modkit_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    echo -e "base\\tmod_code\\tfraction_modified" > ${meta.id}_modkit_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: "stub"
    END_VERSIONS
    """
}

/*
============================================================================================
    MODULE: MODKIT_PILEUP
    Tool   : modkit (>=0.5.0)
    Purpose: Extract per-site methylation pileups from BAM with MM/ML tags.
             Produces bedMethyl format for downstream DMR analysis.
============================================================================================
*/

process MODKIT_PILEUP {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/methylation/${meta.id}", mode: 'copy'

    container   'ontresearch/modkit:latest'

    input:
    tuple val(meta), path(bam), path(bai)
    path  reference
    path  reference_fai
    val   motifs
    val   cpg_only

    output:
    tuple val(meta), path("*.bedmethyl.gz"),   emit: bedmethyl
    tuple val(meta), path("*.bedmethyl.gz.tbi"), emit: tbi
    path  "versions.yml",                      emit: versions

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def cpg_filter  = cpg_only ? '--cpg' : ''
    """
    modkit pileup \\
        --threads ${task.cpus} \\
        ${cpg_filter} \\
        --motif ${motifs} \\
        --reference ${reference} \\
        --output ${prefix}.bedmethyl.gz \\
        ${bam}

    tabix -p vcf ${prefix}.bedmethyl.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version 2>&1 | awk '{print \$2}')
        tabix: \$(tabix --version 2>&1 | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bedmethyl.gz
    touch ${meta.id}.bedmethyl.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: MODKIT_SAMPLE_PROBS
    Purpose: Inspect distribution of modification probability scores per sample.
             Used to assess whether the model threshold (0.5) is appropriate and
             whether 5hmC is separable from 5mC in this dataset.
========================================================================================
*/

process MODKIT_SAMPLE_PROBS {
    tag         "$meta.id"
    label       'process_low'
    publishDir  "${params.outdir}/qc/methylation/${meta.id}", mode: 'copy'

    container   'ontresearch/modkit:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_prob_summary.tsv"), emit: summary
    path  "versions.yml",                                 emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    modkit sample-probs \\
        --threads ${task.cpus} \\
        --only-mapped \\
        --out-dir . \\
        --prefix ${prefix} \\
        ${bam}

    mv ${prefix}*.tsv ${prefix}_prob_summary.tsv 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_prob_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: "stub"
    END_VERSIONS
    """
}

/*
========================================================================================
    MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    Purpose: Aggregate all versions.yml files into a single MultiQC-compatible YAML.
========================================================================================
*/

process CUSTOM_DUMPSOFTWAREVERSIONS {
    label       'process_low'

    container   'community.wave.seqera.io/library/python:3.11'

    input:
    path versions

    output:
    path "software_versions.yml",      emit: yml
    path "software_versions_mqc.yml",  emit: mqc_yml
    path "versions.yml",               emit: versions

    script:
    """
    #!/usr/bin/env python3
    import yaml, glob

    all_versions = {}
    for f in glob.glob('*.yml'):
        if f == 'versions.yml': continue
        with open(f) as fh:
            data = yaml.safe_load(fh) or {}
            all_versions.update(data)

    with open('software_versions.yml', 'w') as out:
        yaml.dump(all_versions, out, default_flow_style=False)

    mqc = {
        'id': 'software_versions',
        'section_name': 'Software Versions',
        'section_href': 'https://github.com/your-org/ont-pipeline',
        'plot_type': 'html',
        'description': 'Tool versions used in this pipeline run.',
        'data': '<dl class="dl-horizontal">' +
                ''.join(f'<dt>{p}</dt><dd><samp>{v}</samp></dd>'
                        for sect in all_versions.values()
                        for p, v in (sect.items() if isinstance(sect, dict) else {}).items()) +
                '</dl>'
    }
    with open('software_versions_mqc.yml', 'w') as out:
        yaml.dump(mqc, out, default_flow_style=False)

    with open('versions.yml', 'w') as out:
        out.write('"${task.process}":\\n  python: stub\\n')
    """

    stub:
    """
    touch software_versions.yml
    touch software_versions_mqc.yml
    touch versions.yml
    """
}

/*
============================================================================================
    MODULE: DORADO_BASECALLER
    Tool   : Dorado (>=0.5.0)
    Purpose: Basecall raw pod5/fast5 signal → FASTQ or uBAM with MM/ML tags.
             Uses SUP model by default for highest accuracy.
============================================================================================
*/

process DORADO_BASECALLER {
    tag         "$meta.id"
    label       'process_gpu'

    container   'community.wave.seqera.io/library/dorado:latest'

    input:
    tuple val(meta), path(pod5)
    val   model
    val   min_qscore

    output:
    tuple val(meta), path("*.ubam"), emit: ubam
    path  "versions.yml",            emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dorado basecaller \\
        ${model} \\
        ${pod5} \\
        --emit-moves \\
        --min-qscore ${min_qscore} \\
        --threads ${task.cpus} \\
        --output-dir . \\
        --ubam \\
        --no-progress

    mv *.ubam ${prefix}.ubam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.ubam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: "stub"
    END_VERSIONS
    """
}

/*
============================================================================================
    MODULE: DORADO_DUPLEX
    Tool   : Dorado (>=0.5.0)
    Purpose: Perform duplex basecalling (template + complement pairing).
             Produces higher accuracy reads (Q50+) by combining forward/reverse strands.
============================================================================================
*/

process DORADO_DUPLEX {
    tag         "$meta.id"
    label       'process_gpu'

    container   'community.wave.seqera.io/library/dorado:latest'

    input:
    tuple val(meta), path(pod5)
    val   model

    output:
    tuple val(meta), path("*.ubam"), emit: ubam
    path  "versions.yml",            emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dorado duplex \\
        ${model} \\
        ${pod5} \\
        --emit-moves \\
        --threads ${task.cpus} \\
        --output-dir . \\
        --ubam \\
        --no-progress

    mv *.ubam ${prefix}.duplex.ubam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.duplex.ubam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: "stub"
    END_VERSIONS
    """
}

/*
============================================================================================
    MODULE: DORADO_MERGE
    Tool   : Dorado (>=0.5.0)
    Purpose: Merge simplex and duplex uBAMs, removing duplicate reads.
             Duplex reads supersede their simplex counterparts.
============================================================================================
*/

process DORADO_MERGE {
    tag         "$meta.id"
    label       'process_low'

    container   'community.wave.seqera.io/library/dorado:latest'

    input:
    tuple val(meta), path(simplex_ubam), path(duplex_ubam)

    output:
    tuple val(meta), path("*.ubam"), emit: ubam
    path  "versions.yml",             emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dorado merge \\
        ${simplex_ubam} \\
        ${duplex_ubam} \\
        --output-dir . \\
        --ubam \\
        --no-progress

    mv *.ubam ${prefix}.merged.ubam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.merged.ubam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: "stub"
    END_VERSIONS
    """
}

/*
============================================================================================
    MODULE: MINIMAP2_ALIGN
    Tool   : minimap2 (>=2.26)
    Purpose: Align long reads (ONT/PacBio) to reference genome.
             Uses lr:hq preset for R10.4.1 SUP reads.
============================================================================================
*/

process MINIMAP2_ALIGN {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/minimap2:latest'

    input:
    tuple val(meta), path(reads)
    path  reference
    path  reference_fai
    val   preset
    val   extra_flags

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml",            emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def is_ubam  = reads.toString().endsWith('.bam')
    def copy_tags = is_ubam ? '-y' : ''
    """
    minimap2 \\
        -t ${task.cpus} \\
        -ax ${preset} \\
        ${copy_tags} \\
        ${extra_flags} \\
        ${reference} \\
        ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.bam -

    samtools index -@ ${task.cpus} ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1 | head -1 | awk '{print \$2}')
        samtools: \$(samtools --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """
}

/*
============================================================================================
    MODULE: SEQKIT_VALIDATE
    Tool   : seqkit (>=2.0)
    Purpose: Validate FASTQ files and remove malformed reads where SEQ/QUAL
             lengths differ. Prevents minimap2/samtools from failing with
             "SEQ and QUAL are of different length" errors.
============================================================================================
*/

process SEQKIT_VALIDATE {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/seqkit:latest'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.validated.fastq.gz"), emit: reads
    path  "versions.yml",                                      emit: versions

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def decompress = reads.toString().endsWith('.gz') ? 'zcat' : 'cat'
    """
    ${decompress} ${reads} \\
        | seqkit seq -v \\
        | gzip > ${prefix}.validated.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.validated.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: "stub"
    END_VERSIONS
    """
}

