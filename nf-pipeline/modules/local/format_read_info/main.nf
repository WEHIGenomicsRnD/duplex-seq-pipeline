process FORMAT_READ_INFO {
    tag "$meta.id"
    label 'process_medium'

    conda 'conda-forge::python=3.12 conda-forge::numpy conda-forge::pandas'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2571c2-0' :
        'biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2571c2-0' }"

    input:
    tuple val(meta), path(rinfo)

    output:
    tuple val(meta), path("*.txt.gz"), emit: formatted
    path 'versions.yml',               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    format_read_info.py ${rinfo} ${prefix}.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
