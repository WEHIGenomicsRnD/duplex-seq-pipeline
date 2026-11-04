process CALC_DUPLEX_METRICS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/wehigenomicsrnd/calculate-duplex-metrics:0.2.0'

    input:
    tuple val(meta), path(input)
    // mapped BAM, or the read info table from GET_READ_INFO
    tuple val(meta2), path(fasta)
    path target_bed
    val input_format

    output:
    tuple val(meta), path("*_metrics.csv"), emit: csv
    tuple val("${task.process}"), val('calc-duplex-metrics'), eval("echo '${task.container}' | sed 's/^.*://'"), topic: versions, emit: versions_calcduplexmetrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // ext.args2 (--sort_mem, --bam_chunk_size) and --cores only apply to BAM input, as in Snakemake
    def args2 = input_format == 'bam' ? "--cores ${task.cpus} ${task.ext.args2 ?: ''}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def target = target_bed ? "--target_bed ${target_bed}" : ''
    """
    calc-duplex-metrics \\
        --input ${input} \\
        --input_format ${input_format} \\
        --output ${prefix}_metrics.csv \\
        --ref_fasta ${fasta} \\
        --sample ${meta.id} \\
        ${target} \\
        ${args2} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample,metric,value\\n' > ${prefix}_metrics.csv
    """
}
