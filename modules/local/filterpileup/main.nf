process FILTER_PILEUP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // filter_pileup.py only needs python; this is the image of FILTER_VCF
    container 'community.wave.seqera.io/library/python_numpy_pandas_pysam:d8e93dd4afe5255c'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_filtered.vcf"), emit: vcf
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // the minimum supporting depth and the ploidy are set in conf/modules.config
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_pileup.py ${args} ${vcf} ${prefix}_filtered.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '##fileformat=VCFv4.3\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}\\n' > ${prefix}_filtered.vcf
    """
}
