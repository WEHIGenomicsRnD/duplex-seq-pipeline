process FGBIO_EXTRACTUMISFROMBAM {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4d/4d1150a2e123f49f8c268f0ab429847afae642376fa52af713b846b084df4a9f/data'
        : 'community.wave.seqera.io/library/fgbio:3.1.2--6e9400d507a9dc55'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('fgbio'), eval('fgbio --version 2>&1 | tr -d "[:cntrl:]" | sed -e "s/^.*Version: //;s/\\[.*$//"'), topic: versions, emit: versions_fgbio

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.mega * 0.9).intValue() : 4096
    """
    fgbio \\
        -Xmx${avail_mem}m \\
        -Djava.io.tmpdir=. \\
        ExtractUmisFromBam \\
        --input ${bam} \\
        --output ${prefix}.bam \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
