process GET_READ_INFO {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/python_numpy_pandas_pysam:d8e93dd4afe5255c'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt.gz"), emit: rinfo
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval('python3 -c "import pysam; print(pysam.__version__)" 2>&1 | tail -1'), topic: versions, emit: versions_pysam
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)" 2>&1 | tail -1'), topic: versions, emit: versions_pandas

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args are extra options for `sort` (the memory buffer, -S, is set in conf/modules.config)
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_read_info.py ${bam} \\
        | sort -k 2,7 -k 1,1 -u --parallel=${task.cpus} ${args} -T . \\
        | cut -f 2-7 \\
        | uniq --count > ${prefix}.rinfo_summary.tmp

    if [ -s ${prefix}.rinfo_summary.tmp ]; then
        format_read_info.py ${prefix}.rinfo_summary.tmp ${prefix}.txt.gz
    else
        echo "No reads found for sample ${meta.id}" >&2
        exit 1
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.txt.gz
    """
}
