rule convert_to_unmapped_bam:
    input:
        r1='fastq/{sample}_R1.fastq.gz',
        r2='fastq/{sample}_R2.fastq.gz'
    output:
        'results/unmapped_bams/{sample}_unmapped.bam',
    log:
        'logs/convert_to_unmapped_bam_{sample}.log'
    envmodules:
        'picard-tools/2.20.2'
    conda:
        '../envs/picard.yaml'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='0-12:0:0'
    shell:
        '''
        FastqToSam \
            F1={input.r1} \
            F2={input.r2} \
            O={output} \
            SM={wildcards.sample} \
            RG=null
        '''
