process FILTER_VCF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // filter_vcf.py only needs python; this is the image of GET_READ_INFO
    container 'community.wave.seqera.io/library/python_numpy_pandas_pysam:d8e93dd4afe5255c'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_filtered.vcf"), emit: vcf
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_vcf.py ${vcf} ${prefix}_filtered.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp ${vcf} ${prefix}_filtered.vcf
    """
}
