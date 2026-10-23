process FILTER_VCF {
    tag "$meta.id"
    label 'process_low'

    conda 'conda-forge::python=3.12'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_filtered.vcf"), emit: vcf
    path 'versions.yml',                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_vcf.py ${vcf} ${prefix}_filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
