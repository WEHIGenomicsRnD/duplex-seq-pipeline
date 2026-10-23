process CALC_DUPLEX_METRICS {
    tag "$meta.id"
    label 'process_medium'

    conda 'conda-forge::r-base=4.4.1 conda-forge::r-remotes conda-forge::r-devtools conda-forge::r-argparse conda-forge::r-data.table conda-forge::r-magrittr bioconda::bioconductor-genomeinfodb bioconda::bioconductor-rsamtools bioconda::bioconductor-biostrings bioconda::bioconductor-genomicranges bioconda::bioconductor-iranges conda-forge::r-seqinr'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-XXX:latest' :
        'biocontainers/mulled-v2-XXX:latest' }"

    input:
    tuple val(meta), path(rinfo)
    path ref

    output:
    tuple val(meta), path("*_metrics.csv"), emit: metrics
    path 'versions.yml',                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metrics = params.duplex_metrics ?: 'all'
    """
    # Install CalcDuplexMetrics if not present
    if ! Rscript --vanilla -e "library(CalcDuplexMetrics)" 2>/dev/null; then
        Rscript --vanilla \\
            -e "devtools::install_github('WEHIGenomicsRnD/calculate-duplex-metrics')"
        SCRIPT_PATH=\$(Rscript --vanilla \\
            -e "cat(system.file('exec', 'calc-duplex-metrics', package='CalcDuplexMetrics'))" \\
            2>/dev/null)
        ln -sf "\$SCRIPT_PATH" "\$CONDA_PREFIX/bin/calc-duplex-metrics"
    fi

    calc-duplex-metrics \\
        --input ${rinfo} \\
        --output ${prefix}_metrics.csv \\
        --ref_fasta ${ref} \\
        --metrics ${metrics}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/.*version //;s/ .*//')
    END_VERSIONS
    """
}
