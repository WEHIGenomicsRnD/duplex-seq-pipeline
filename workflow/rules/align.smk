rule sam_to_fastq:
    input:
        'results/unmapped_bams/{sample}_unmapped_wUMI.bam'
    output:
        'results/unmapped_bams/{sample}.fastq'
    log:
        'logs/sam_to_fastq_{sample}.log'
    conda:
        '../envs/picard.yaml'
    envmodules:
        'picard-tools/2.20.2'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='0-12:0:0'
    shell:
        '''
        SamToFastq \
            I={input} \
            F={output} \
            INTERLEAVE=true
        '''

rule align:
    input:
        'results/unmapped_bams/{sample}.fastq'
    output:
        'results/mapped_bams/{sample}_mapped.bam'
    log:
        'logs/align_{sample}.log'
    conda:
        '../envs/bwa_samtools.yaml'
    envmodules:
        'bwa/0.7.17'
    threads:
        16
    resources:
        mem_mb=36864,
        runtime='0-12:0:0'
    shell:
        '''
        bwa mem -p -t {threads} {config[ref]} {input} | samtools view -bS - > {output}
        '''

rule merge_bam_alignment:
    input:
        unmapped='results/unmapped_bams/{sample}_unmapped_wUMI.bam',
        mapped='results/mapped_bams/{sample}_mapped.bam'
    output:
        'results/mapped_bams/{sample}_mapped_merged.bam'
    log:
        'logs/merge_bam_alignment_{sample}.log'
    conda:
        '../envs/picard.yaml'
    envmodules:
        'picard-tools/2.20.2'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='0-12:0:0'
    shell:
        '''
        MergeBamAlignment \
                UNMAPPED={input.unmapped} \
                ALIGNED={input.mapped} \
                O={output} \
                R={config[ref]} \
                SO=coordinate \
                ALIGNER_PROPER_PAIR_FLAGS=true \
                MAX_GAPS=-1 \
                ORIENTATIONS=FR \
                CREATE_INDEX=true
        '''

rule remap_sam_to_fastq:
    input:
        'results/consensus/{sample}_consensus_unmapped.bam'
    output:
        'results/consensus/{sample}_consensus.fastq'
    log:
        'logs/remap_sam_to_fastq_{sample}.log'
    conda:
        '../envs/picard.yaml'
    envmodules:
        'picard-tools/2.20.2'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='0-12:0:0'
    shell:
        '''
        SamToFastq \
            I={input} \
            F={output} \
            INTERLEAVE=true
        '''

rule remap_align:
    input:
        'results/consensus/{sample}_consensus.fastq'
    output:
        'results/consensus/{sample}_consensus_mapped.bam'
    log:
        'logs/remap_align_{sample}.log'
    conda:
        '../envs/bwa_samtools.yaml'
    envmodules:
        'bwa/0.7.17'
    threads:
        16
    resources:
        mem_mb=36864,
        runtime='0-12:0:0'
    shell:
        '''
        bwa mem -p -t {threads} {config[ref]} {input} | samtools view -bS - > {output}
        '''

rule remap_merge_bam_alignment:
    input:
        unmapped='results/consensus/{sample}_consensus_unmapped.bam',
        mapped='results/consensus/{sample}_consensus_mapped.bam'
    output:
        'results/consensus/{sample}_consensus_mapped_merged.bam'
    log:
        'logs/remap_merge_bam_alignment_{sample}.log'
    conda:
        '../envs/picard.yaml'
    envmodules:
        'picard-tools/2.20.2'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='0-12:0:0'
    shell:
        '''
        MergeBamAlignment \
                UNMAPPED={input.unmapped} \
                ALIGNED={input.mapped} \
                O={output} \
                R={config[ref]} \
                SO=coordinate \
                ALIGNER_PROPER_PAIR_FLAGS=true \
                MAX_GAPS=-1 \
                ORIENTATIONS=FR \
                CREATE_INDEX=true
        '''

