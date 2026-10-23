process FASTQ_SCREEN {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::fastq-screen=0.16.0 bioconda::bowtie2=2.5.4'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20571b:13a668ea63c3900b6880c1539ed56a78abccd960-0' :
        'biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20571b:13a668ea63c3900b6880c1539ed56a78abccd960-0' }"

    input:
    tuple val(meta), path(reads)
    path config_file

    output:
    tuple val(meta), path("*_screen.html"), emit: html
    tuple val(meta), path("*_screen.txt"),  emit: txt
    path 'versions.yml',                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fastq_screen \\
        --conf ${config_file} \\
        --outdir . \\
        --force \\
        --threads ${task.cpus} \\
        ${args} \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq_screen: \$(fastq_screen --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+')
    END_VERSIONS
    """
}
