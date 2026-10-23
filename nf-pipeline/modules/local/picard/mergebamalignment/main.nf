process PICARD_MERGEBAMALIGNMENT {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::picard=3.4.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.4.0--hdfd78af_0' :
        'biocontainers/picard:3.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(unmapped_bam), path(mapped_bam)
    path ref
    path dict

    output:
    tuple val(meta), path("*_merged.bam"), emit: bam
    path 'versions.yml',                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = (task.memory.toMega() * 0.9).intValue()
    """
    picard -Xmx${avail_mem}m MergeBamAlignment \\
        UNMAPPED=${unmapped_bam} \\
        ALIGNED=${mapped_bam} \\
        O=${prefix}_merged.bam \\
        R=${ref} \\
        SO=coordinate \\
        ALIGNER_PROPER_PAIR_FLAGS=true \\
        MAX_GAPS=-1 \\
        ORIENTATIONS=FR \\
        CREATE_INDEX=true \\
        TMP_DIR=\${TMPDIR:-/tmp} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CreateSequenceDictionary --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
