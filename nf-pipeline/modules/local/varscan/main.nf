process VARSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda 'bioconda::samtools=1.21 bioconda::varscan=2.4.6'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-67e53cc0e67f27a56b0eb0d49d51e9c37e223a88:latest' :
        'biocontainers/mulled-v2-67e53cc0e67f27a56b0eb0d49d51e9c37e223a88:latest' }"

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
    def min_coverage = params.min_coverage
    def min_reads2 = params.min_reads2
    def min_vaf = params.min_vaf
    """
    samtools mpileup -f ${ref} ${bam} | \\
        varscan mpileup2snp \\
        --min-avg-qual 0 \\
        --min-coverage ${min_coverage} \\
        --min-reads2 ${min_reads2} \\
        --min-var-freq ${min_vaf} \\
        --strand-filter 0 \\
        --p-value 1 \\
        --output-vcf 1 > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | head -n1 | sed 's/samtools //')
        varscan: \$(varscan 2>&1 | grep -oP 'VarScan v\\K[\\d.]+')
    END_VERSIONS
    """
}
