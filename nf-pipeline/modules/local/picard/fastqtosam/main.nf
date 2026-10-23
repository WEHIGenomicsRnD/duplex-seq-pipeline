process PICARD_FASTQTOSAM {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::picard=3.4.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.4.0--hdfd78af_0' :
        'biocontainers/picard:3.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_unmapped.bam"), emit: bam
    path 'versions.yml',                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def (r1, r2) = reads instanceof List ? reads : [reads, null]
    def fastq2_arg = r2 ? "F2=${r2}" : ''
    def avail_mem = (task.memory.toMega() * 0.9).intValue()
    """
    picard -Xmx${avail_mem}m FastqToSam \\
        F1=${r1} \\
        ${fastq2_arg} \\
        O=${prefix}_unmapped.bam \\
        SM=${meta.id} \\
        RG=null \\
        TMP_DIR=\${TMPDIR:-/tmp} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CreateSequenceDictionary --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
