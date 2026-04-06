/*
=======================================================================================
    MODULE: BGZIP_TABIX
    Tool   : htslib (>=1.18)
    Purpose: Compress VCF/BED files with bgzip and index with tabix.
=======================================================================================
*/

process BGZIP_TABIX {
    tag         "$meta.id"
    label       'process_low'

    container   'sktrinh12/samtools:latest'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("${meta.id}.out.gz"), emit: compressed
    tuple val(meta), path("${meta.id}.out.gz.tbi"), emit: tbi, optional: true
    path  "versions.yml",                                         emit: versions

    script:
    def is_compressed = file.toString() =~ /\.vcf(\.gz)?$/ || file.toString() =~ /\.bed(\.gz)?$/
    def preset = file.toString() =~ /\.vcf(\.gz)?$/ ? 'vcf' : 
                 file.toString() =~ /\.bed(\.gz)?$/ ? 'bed' : 'vcf'
    """
    realpath_cmd=\$(realpath ${file})
    if [ "${is_compressed}" == "true" ]; then
        # Already compressed, symlink to new name and index
        ln -sf "\$realpath_cmd" ${meta.id}.out.gz
        /opt/conda/bin/tabix --force -p ${preset} ${meta.id}.out.gz
    else
        # Compress and index
        /opt/conda/bin/bgzip -@ ${task.cpus} -c ${file} > ${meta.id}.out.gz
        /opt/conda/bin/tabix -p ${preset} ${meta.id}.out.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "BGZIP_TABIX":
        htslib: \$(/opt/conda/bin/bgzip --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.out.gz
    touch ${meta.id}.out.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "BGZIP_TABIX":
        htslib: "stub"
    END_VERSIONS
    """
}
