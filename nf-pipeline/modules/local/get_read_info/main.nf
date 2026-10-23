process GET_READ_INFO {
    tag "$meta.id"
    label 'process_high'

    conda 'conda-forge::python=3.12 conda-forge::numpy conda-forge::pandas bioconda::pysam=0.22.1'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2fe2a8bea2d45e4a85a199e3aaa700d4ee9bf164:af4d7dccb01cf5bb5c1b5a4cf59b97a3f3b9debc-0' :
        'biocontainers/mulled-v2-2fe2a8bea2d45e4a85a199e3aaa700d4ee9bf164:af4d7dccb01cf5bb5c1b5a4cf59b97a3f3b9debc-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_rinfo.txt"), emit: rinfo
    path 'versions.yml',                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def buffer_size = task.ext.buffer_size ?: '50G'
    """
    get_read_info.py ${bam} | \\
        sort -k 2,6 -k 1,1 -u --parallel=${task.cpus} -S ${buffer_size} -T \${TMPDIR:-/tmp} | \\
        cut -f 2-6 | \\
        uniq --count > ${prefix}_rinfo.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        pysam: \$(python -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """
}
