process CALC_DUPLEX_METRICS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/wehigenomicsrnd/calculate-duplex-metrics:0.2.0'

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    path target_bed

    output:
    tuple val(meta), path("*_metrics.csv"), emit: csv
    tuple val("${task.process}"), val('calc-duplex-metrics'), eval("echo '${task.container}' | sed 's/^.*://'"), topic: versions, emit: versions_calcduplexmetrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def target = target_bed ? "--target_bed ${target_bed}" : ''
    """
    calc-duplex-metrics \\
        --input ${bam} \\
        --input_format bam \\
        --output ${prefix}_metrics.csv \\
        --ref_fasta ${fasta} \\
        --sample ${meta.id} \\
        --cores ${task.cpus} \\
        ${target} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample,metric,value\\n' > ${prefix}_metrics.csv
    """
}
