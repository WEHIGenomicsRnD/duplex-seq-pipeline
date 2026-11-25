process NVC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/nvc:0.0.4--cd3232ac269fda45'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('nvc'), eval('naive_variant_caller.py --version 2>&1 | sed "s/^.*version //"'), topic: versions, emit: versions_nvc

    when:
    task.ext.when == null || task.ext.when

    script:
    // the caller's options (--ploidy, --min_support_depth, ...) are set in conf/modules.config
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    naive_variant_caller.py \\
        --bam=${bam} \\
        --output_vcf_filename=${prefix}.vcf \\
        --reference_genome_filename=${fasta} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' > ${prefix}.vcf
    """
}
