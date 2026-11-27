process BCFTOOLS_PILEUP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // the bcftools image of the nf-core bcftools/mpileup module
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/914b9b38bbfa22c458ceee6c578fd3d70ee1bffc33d3a064b0663ed4944e7872/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.25--130aa7603bf3bf5c'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)
    path regions

    output:
    tuple val(meta), path("*.pileup.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    // args: the mpileup options; args2: the bcftools view filter that selects the candidate positions
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_cmd = regions ? "--regions-file ${regions}" : ''
    """
    # the pileup of every covered position, from which all records (nucleotides and indels) of
    # the candidate positions are kept
    bcftools mpileup \\
        --fasta-ref ${fasta} \\
        ${regions_cmd} \\
        --annotate FORMAT/AD,FORMAT/DP \\
        ${args} \\
        --output-type b \\
        --output ${prefix}.all.bcf \\
        ${bam}

    bcftools view ${args2} ${prefix}.all.bcf \\
        | bcftools query --format '%CHROM\\t%POS\\n' \\
        > candidate_positions.txt

    if [ -s candidate_positions.txt ]; then
        bcftools view --targets-file candidate_positions.txt --output-type z --output ${prefix}.pileup.vcf.gz ${prefix}.all.bcf
    else
        bcftools view --header-only --output-type z --output ${prefix}.pileup.vcf.gz ${prefix}.all.bcf
    fi
    rm ${prefix}.all.bcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' | gzip > ${prefix}.pileup.vcf.gz
    """
}
