process NVC {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::nvc=0.0.4 bioconda::samtools=1.21'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-XXX:latest' :
        'biocontainers/mulled-v2-XXX:latest' }"

    input:
    tuple val(meta), path(bam)
    path ref

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path 'versions.yml',            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ploidy = params.ploidy
    def min_reads2 = params.min_reads2
    def min_mapq = params.min_mapq
    """
    samtools index ${bam}

    naive_variant_caller.py \\
        --bam=${bam} \\
        --output_vcf_filename=${prefix}.vcf \\
        --reference_genome_filename=${ref} \\
        --variants_only \\
        --ploidy=${ploidy} \\
        --min_support_depth=${min_reads2} \\
        --min_mapping_quality=${min_mapq} \\
        --min_base_quality=0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nvc: \$(naive_variant_caller.py --version 2>&1 || echo '0.0.4')
        samtools: \$(samtools version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
