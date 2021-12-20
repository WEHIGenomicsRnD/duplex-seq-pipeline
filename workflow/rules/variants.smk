rule clip_bam:
    input:
        'results/consensus/{sample}_consensus_mapped_merged_filtered.bam'
    output:
        'results/consensus/{sample}_consensus_mapped_merged_filtered_clipped.bam'
    log:
        'logs/clip_bam_{sample}.log'
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
            ClipBam \
            --input={input} \
            --output={output} \
            --ref {config[ref]} \
            --clipping-mode=Hard \
            --clip-overlapping-reads=true
        '''

rule varscan:
    input:
        'results/consensus/{sample}_consensus_mapped_merged_filtered_clipped.bam'
    output:
        'results/variants/{sample}.vcf'
    log:
        'logs/varscan_{sample}.log'
    envmodules:
        'samtools/1.9',
        'varscan/2.3.9'
    conda:
        '../envs/varscan.yaml'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='1-0:0:0'
    shell:
        '''
        samtools mpileup -f {config[ref]} {input} | \
            varscan mpileup2snp -Xmx{resources.mem_mb} \
            --min-coverage 1 \
            --min-reads2 1 \
            --min-var-freq {config[min_vaf]} \
            --p-value {config[pvalue]} \
            --output-vcf 1 > {output}
        '''
