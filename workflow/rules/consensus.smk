rule call_duplex_consensus:
    input:
        'results/grouped_by_umi/{sample}_grouped.bam'
    output:
        'results/consensus/{sample}_consensus_unmapped.bam'
    log:
        'logs/call_duplex_consensus_{sample}.log'
    conda:
        '../envs/fgbio.yaml'
    threads:
        1
    resources:
        mem_mb=65536,
        runtime='0-12:0:0'
    shell:
        '''
        fgbio -Xmx{resources.mem_mb}m \
            CallDuplexConsensusReads \
            --input={input} \
            --output={output} \
            --error-rate-pre-umi={config[error_rate_pre_umi]} \
            --error-rate-post-umi={config[error_rate_post_umi]} \
            --min-input-base-quality={config[min_input_baseq]}
        '''

rule filter_consensus_reads:
    input:
        'results/consensus/{sample}_consensus_mapped_merged.bam'
    output:
        'results/consensus/{sample}_consensus_mapped_merged_filtered.bam'
    log:
        'logs/filter_consensus_reads_{sample}.log'
    conda:
        '../envs/fgbio.yaml'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='0-12:0:0'
    shell:
        '''
        fgbio -Xmx{resources.mem_mb}m \
            FilterConsensusReads \
            --input={input} \
            --output={output} \
            --ref {config[ref]} \
            --min-reads={config[min_reads]} \
            --max-read-error-rate={config[max_read_error_rate]} \
            --max-base-error-rate={config[max_base_error_rate]} \
            --min-base-quality={config[min_base_quality]} \
            --max-no-call-fraction={config[max_no_call_fraction]}
        '''
