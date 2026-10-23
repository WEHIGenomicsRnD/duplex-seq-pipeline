process FGBIO_EXTRACT_UMIS {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::fgbio=2.5.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.5.0--hdfd78af_0' :
        'biocontainers/fgbio:2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_wUMI.bam"), emit: bam
    path 'versions.yml',                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_structure = params.read_structure
    def mol_tags = params.molecular_index_tags
    def single_tag = params.single_tag
    """
    fgbio -Xmx${task.memory.toMega()}m \\
        -Djava.io.tmpdir=\${TMPDIR:-/tmp} \\
        ExtractUmisFromBam \\
        --input=${bam} \\
        --output=${prefix}_wUMI.bam \\
        --read-structure=${read_structure} \\
        --molecular-index-tags=${mol_tags} \\
        --single-tag=${single_tag} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | head -n1)
    END_VERSIONS
    """
}
