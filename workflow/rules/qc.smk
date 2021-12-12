rule fastQC:
    input:
        r1='fastq/{sample}_R1.fastq.gz',
        r2='fastq/{sample}_R2.fastq.gz',
    conda:
        '../envs/fastqc.yaml'
    envmodules:
        'fastqc/0.11.8',
    threads:
        2
    resources:
        mem_mb=8192,
        runtime='0-12:0:0',
    log:
        'logs/fastQC_{sample}.log'
    output:
        report1='results/QC/fastQC/{sample}_R1_fastqc.html',
        report2='results/QC/fastQC/{sample}_R2_fastqc.html',
    shell:
        '''
        fastqc -t {threads} {input.r1} {input.r2} -o results/QC/fastQC
        '''

rule fastqScreen:
    input:
        r1='fastq/{sample}_R1.fastq.gz',
        r2='fastq/{sample}_R2.fastq.gz',
        config_file=config['FastqScreen_config'],
    conda:
        '../envs/fastqscreen.yaml'
    envmodules:
        'bowtie2/2.3.4.1',
    threads:
        8
    resources:
        mem_mb=8192,
        runtime='1-0:0:0',
    log:
        'logs/fastqScreen_{sample}.log'
    output:
        report3='results/QC/fastqScreen/{sample}_R1_screen.html',
        report4='results/QC/fastqScreen/{sample}_R2_screen.html'
    shell:
        '''
        fastq_screen --conf {input.config_file} \
            --outdir results/QC/fastqScreen \
            {input.r1} {input.r2}
        '''

rule samtools_stats:
    input:
        'results/mapped_bams/{sample}_mapped_merged.bam'
    output:
        stats='results/QC/samtools/stats/{sample}.txt',
        flagstat='results/QC/samtools/flagstat/{sample}.txt',
    log:
        'logs/samtools_stats_{sample}.log'
    envmodules:
        'samtools/1.19'
    conda:
        '../envs/bwa_samtools.yaml'
    threads:
        1
    resources:
        mem_mb=8192,
        runtime='0-5:0:0'
    shell:
        '''
        samtools stats {input} > {output.stats}
        samtools flagstat {input} > {output.flagstat}
        '''

rule qualimap:
    input:
        'results/mapped_bams/{sample}_mapped_merged.bam'
    output:
        'results/QC/qualimap/{sample}/qualimapReport.html'
    log:
        'logs/qualimap_{sample}.log'
    conda:
        '../envs/qualimap.yaml'
    threads:
        8
    resources:
        mem_mb=8192,
        runtime='0-6:0:0'
    shell:
        '''
        qualimap bamqc \
            --java-mem-size={resources.mem_mb}M \
            -bam {input} \
            -outdir results/QC/qualimap/{wildcards.sample} \
            -outformat HTML \
            -nt {threads}
        '''

rule multiQC:
    input:
        fastqc=expand(
            'results/QC/fastQC/{sample}_R1_fastqc.html',
            sample=samples
        ),
        fastqscreen=expand(
            'results/QC/fastqScreen/{sample}_R1_screen.html',
            sample=samples
        ),
        qualimap=expand(
            'results/QC/qualimap/{sample}/qualimapReport.html',
            sample=samples
        ),
        samtools=expand(
            'results/QC/samtools/stats/{sample}.txt',
            sample=samples
        )
    conda:
        '../envs/multiqc.yaml'
    envmodules:
        'MultiQC/1.10.1',
    resources:
        mem_mb=65536,
        runtime='0-12:0:0'
    log:
        'logs/multiQC.log'
    output:
        'results/QC/multiQC/multiqc_report.html'
    shell:
        '''
        multiqc results/QC/fastQC/ \
            results/QC/fastqScreen/ \
            results/QC/qualimap/ \
            results/QC/samtools/ \
            -o results/QC/multiQC -f
        '''

rule samtools_stats_consensus:
    input:
        'results/consensus/{sample}_consensus_mapped_merged_filtered.bam'
    output:
        stats='results/QC/consensus/samtools/stats/{sample}.txt',
        flagstat='results/QC/consensus/samtools/flagstat/{sample}.txt',
    log:
        'logs/samtools_stats_{sample}.log'
    envmodules:
        'samtools/1.19'
    conda:
        '../envs/bwa_samtools.yaml'
    threads:
        1
    resources:
        mem_mb=8192,
        runtime='0-5:0:0'
    shell:
        '''
        samtools stats {input} > {output.stats}
        samtools flagstat {input} > {output.flagstat}
        '''

rule qualimap_consensus:
    input:
        'results/consensus/{sample}_consensus_mapped_merged_filtered.bam'
    output:
        'results/QC/consensus/qualimap/{sample}/qualimapReport.html'
    log:
        'logs/qualimap_consensus_{sample}.log'
    conda:
        '../envs/qualimap.yaml'
    threads:
        8
    resources:
        mem_mb=8192,
        runtime='0-6:0:0'
    shell:
        '''
        qualimap bamqc \
            --java-mem-size={resources.mem_mb}M \
            -bam {input} \
            -outdir results/QC/consensus/qualimap/{wildcards.sample} \
            -outformat HTML \
            -nt {threads}
        '''

rule multiQC_consensus:
    input:
        qualimap=expand(
            'results/QC/consensus/qualimap/{sample}/qualimapReport.html',
            sample=samples
        ),
        samtools=expand(
            'results/QC/consensus/samtools/stats/{sample}.txt',
            sample=samples
        )
    conda:
        '../envs/multiqc.yaml'
    envmodules:
        'MultiQC/1.10.1',
    resources:
        mem_mb=65536,
        runtime='0-12:0:0'
    log:
        'logs/multiQC_consensus.log'
    output:
        'results/QC/consensus/multiQC/multiqc_report.html'
    shell:
        '''
        multiqc \
            results/QC/consensus/qualimap/ \
            results/QC/consensus/samtools/ \
            -o results/QC/consensus/multiQC -f
        '''
