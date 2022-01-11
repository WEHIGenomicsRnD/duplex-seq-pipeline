rule assemble:
    input:
        'results/consensus/{sample}_consensus.fastq'
    output:
        'results/QC/checkm/assembly/{sample}.contig'
    log:
        'logs/assemble_{sample}.log'
    conda:
        '../envs/assemble.yaml'
    threads:
        8
    resources:
        mem_mb=65536,
        runtime='1-0:0:0'
    shell:
        '''
        echo -e \
            "max_rd_len={config[max_rd_len]}\n[LIB]\nq={input}" \
            > results/QC/checkm/assembly/{wildcards.sample}.config ;
        SOAPdenovo-127mer all \
            -s results/QC/checkm/assembly/{wildcards.sample}.config \
            -p {threads} \
            -o results/QC/checkm/assembly/{wildcards.sample}
        '''

rule checkm:
    input:
        ['results/QC/checkm/assembly/{sample}.contig'.format(
            sample=sample
         ) for sample in samples]
    output:
        'results/QC/checkm/stats.txt'
    log:
        'logs/checkm.log'
    conda:
        '../envs/checkm.yaml'
    threads:
        8
    resources:
        mem_mb=32768,
        runtime='0-12:0:0'
    shell:
        '''
        checkm lineage_wf \
            -t {threads} \
            -x contig \
            results/QC/checkm/assembly \
            checkm > {output}
        '''
