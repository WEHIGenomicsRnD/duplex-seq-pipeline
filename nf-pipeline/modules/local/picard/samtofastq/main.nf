process PICARD_SAMTOFASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::picard=3.4.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.4.0--hdfd78af_0' :
        'biocontainers/picard:3.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq"), emit: fastq
    path 'versions.yml',              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = (task.memory.toMega() * 0.9).intValue()
    """
    picard -Xmx${avail_mem}m SamToFastq \\
        I=${bam} \\
        F=${prefix}.fastq \\
        INTERLEAVE=true \\
        TMP_DIR=\${TMPDIR:-/tmp} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CreateSequenceDictionary --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
