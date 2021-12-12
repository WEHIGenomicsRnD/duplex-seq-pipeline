rule extract_umis:
    input:
        'results/unmapped_bams/{sample}_unmapped.bam'
    output:
        'results/unmapped_bams/{sample}_unmapped_uncorrected_UMI.bam'
    log:
        'logs/extract_umis_{sample}.log'
    conda:
        '../envs/fgbio.yaml'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='0-12:0:0'
    shell:
        '''
        fgbio ExtractUmisFromBam \
            --input={input} \
            --output={output} \
            --read-structure={config[read_structure]} \
            --molecular-index-tags={config[molecular_index_tags]} \
            --single-tag={config[single_tag]}
        '''

rule correct_umis:
    input:
        bam='results/unmapped_bams/{sample}_unmapped_uncorrected_UMI.bam',
        umis=config['umis']
    output:
        'results/unmapped_bams/{sample}_unmapped_wUMI.bam'
    log:
        'logs/correct_umis_{sample}.log'
    conda:
        '../envs/fgbio.yaml'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='0-12:0:0'
    shell:
        '''
        fgbio CorrectUmis \
            --input={input.bam} \
            --output={output} \
            --max-mismatches={config[max_mismatches]} \
            --min-distance={config[min_distance]} \
            -t RX -U {input.umis}
        '''

rule group_reads_by_umi:
    input:
        'results/mapped_bams/{sample}_mapped_merged.bam'
    output:
        'results/grouped_by_umi/{sample}_grouped.bam'
    log:
        'logs/group_reads_by_umi_{sample}.log'
    conda:
        '../envs/fgbio.yaml'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='0-12:0:0'
    shell:
        '''
        fgbio GroupReadsByUmi \
            --input={input} \
            --output={output} \
            --strategy=paired \
            --edits=1 \
            --min-map-q={config[min_mapq]} \
        '''
