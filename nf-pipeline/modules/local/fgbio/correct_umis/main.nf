process FGBIO_CORRECT_UMIS {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::fgbio=2.5.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.5.0--hdfd78af_0' :
        'biocontainers/fgbio:2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path umis

    output:
    tuple val(meta), path("*_corrected.bam"), emit: bam
    path 'versions.yml',                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def max_mismatches = params.max_mismatches
    def min_distance = params.min_distance
    """
    fgbio -Xmx${task.memory.toMega()}m \\
        -Djava.io.tmpdir=\${TMPDIR:-/tmp} \\
        CorrectUmis \\
        --input=${bam} \\
        --output=${prefix}_corrected.bam \\
        --max-mismatches=${max_mismatches} \\
        --min-distance=${min_distance} \\
        -t RX \\
        -U ${umis} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | head -n1)
    END_VERSIONS
    """
}
