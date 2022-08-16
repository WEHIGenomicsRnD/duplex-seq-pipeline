rule fastQC:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
    output:
        report1="results/QC/fastQC/{sample}_R1_fastqc.html",
        report2="results/QC/fastQC/{sample}_R2_fastqc.html",
    log:
        "logs/fastQC_{sample}.log",
    conda:
        "../envs/fastqc.yaml"
    envmodules:
        "fastqc/0.11.8",
    threads: cluster["fastqc"]["threads"]
    resources:
        mem_mb=cluster["fastqc"]["mem_mb"],
        runtime=cluster["fastqc"]["runtime"],
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} -o results/QC/fastQC
        """


rule fastqScreen:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
        config_file=config["FastqScreen_config"],
    output:
        report3="results/QC/fastqScreen/{sample}_R1_screen.html",
        report4="results/QC/fastqScreen/{sample}_R2_screen.html",
    log:
        "logs/fastqScreen_{sample}.log",
    conda:
        "../envs/fastqscreen.yaml"
    envmodules:
        "bowtie2/2.3.4.1",
    threads: cluster["fastqscreen"]["threads"]
    resources:
        mem_mb=cluster["fastqscreen"]["mem_mb"],
        runtime=cluster["fastqscreen"]["runtime"],
    shell:
        """
        fastq_screen --conf {input.config_file} \
            --outdir results/QC/fastqScreen \
            --force \
            {input.r1} {input.r2}
        """


rule mark_duplicates:
    input:
        "results/mapped_bams/{sample}_mapped_merged.bam",
    output:
        bam="results/mapped_bams/{sample}_mapped_merged_markdup.bam",
        metrics="results/QC/mark_duplicates/{sample}.txt",
    log:
        "logs/mark_duplicates_{sample}.log",
    conda:
        "../envs/picard.yaml"
    threads: cluster["picard"]["threads"]
    resources:
        mem_mb=cluster["picard"]["mem_mb"],
        runtime=cluster["picard"]["runtime"],
    params:
        max_records_in_ram=cluster["picard"]["max_records_in_ram"],
        java_mem=cluster["picard"]["java_mem"],
    shell:
        """
        picard -Xmx{params.java_mem}m MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            BARCODE_TAG=RX \
            TMP_DIR={resources.tmpdir} \
            MAX_RECORDS_IN_RAM={params.max_records_in_ram}
        """


rule samtools_stats:
    input:
        "results/mapped_bams/{sample}_mapped_merged.bam",
    output:
        stats="results/QC/samtools/stats/{sample}.txt",
        flagstat="results/QC/samtools/flagstat/{sample}.txt",
    log:
        "logs/samtools_stats_{sample}.log",
    envmodules:
        "samtools/1.19",
    conda:
        "../envs/bwa_samtools.yaml"
    threads: cluster["samtools_stats"]["threads"]
    resources:
        mem_mb=cluster["samtools_stats"]["mem_mb"],
        runtime=cluster["samtools_stats"]["runtime"],
    shell:
        """
        samtools stats {input} > {output.stats}
        samtools flagstat {input} > {output.flagstat}
        """


rule qualimap:
    input:
        "results/mapped_bams/{sample}_mapped_merged.bam",
    output:
        "results/QC/qualimap/{sample}/qualimapReport.html",
    log:
        "logs/qualimap_{sample}.log",
    conda:
        "../envs/qualimap.yaml"
    threads: cluster["qualimap"]["threads"]
    resources:
        mem_mb=cluster["qualimap"]["mem_mb"],
        runtime=cluster["qualimap"]["runtime"],
    shell:
        """
        qualimap bamqc \
            --java-mem-size={resources.mem_mb}M \
            -bam {input} \
            -outdir results/QC/qualimap/{wildcards.sample} \
            -outformat HTML \
            -nt {threads}
        """


rule multiQC:
    input:
        fastqc=expand("results/QC/fastQC/{sample}_R1_fastqc.html", sample=samples),
        fastqscreen=expand(
            "results/QC/fastqScreen/{sample}_R1_screen.html", sample=samples
        ),
        qualimap=expand(
            "results/QC/qualimap/{sample}/qualimapReport.html", sample=samples
        ),
        samtools=expand("results/QC/samtools/stats/{sample}.txt", sample=samples),
        markdups=expand("results/QC/mark_duplicates/{sample}.txt", sample=samples),
        config_file=config["multiqc_config"],
    output:
        "results/QC/multiQC/multiqc_report.html",
    log:
        "logs/multiQC.log",
    conda:
        "../envs/multiqc.yaml"
    envmodules:
        "MultiQC/1.10.1",
    threads: cluster["multiqc"]["threads"]
    resources:
        mem_mb=cluster["multiqc"]["mem_mb"],
        runtime=cluster["multiqc"]["runtime"],
    shell:
        """
        multiqc \
            --config {input.config_file} \
            results/QC/fastQC/ \
            results/QC/fastqScreen/ \
            results/QC/qualimap/ \
            results/QC/samtools/ \
            results/QC/mark_duplicates \
            -o results/QC/multiQC -f
        """


rule samtools_stats_consensus:
    input:
        "results/consensus/{sample}_consensus_mapped_merged_filtered.bam",
    output:
        stats="results/QC/consensus/samtools/stats/{sample}.txt",
        flagstat="results/QC/consensus/samtools/flagstat/{sample}.txt",
    log:
        "logs/samtools_stats_{sample}.log",
    envmodules:
        "samtools/1.19",
    conda:
        "../envs/bwa_samtools.yaml"
    threads: cluster["samtools_stats"]["threads"]
    resources:
        mem_mb=cluster["samtools_stats"]["mem_mb"],
        runtime=cluster["samtools_stats"]["runtime"],
    shell:
        """
        samtools stats {input} > {output.stats}
        samtools flagstat {input} > {output.flagstat}
        """


rule qualimap_consensus:
    input:
        "results/consensus/{sample}_consensus_mapped_merged_filtered.bam",
    output:
        "results/QC/consensus/qualimap/{sample}/qualimapReport.html",
    log:
        "logs/qualimap_consensus_{sample}.log",
    conda:
        "../envs/qualimap.yaml"
    threads: cluster["qualimap"]["threads"]
    resources:
        mem_mb=cluster["qualimap"]["mem_mb"],
        runtime=cluster["qualimap"]["runtime"],
    shell:
        """
        qualimap bamqc \
            --java-mem-size={resources.mem_mb}M \
            -bam {input} \
            -outdir results/QC/consensus/qualimap/{wildcards.sample} \
            -outformat HTML \
            -nt {threads}
        """


rule multiQC_consensus:
    input:
        qualimap=expand(
            "results/QC/consensus/qualimap/{sample}/qualimapReport.html",
            sample=samples,
        ),
        samtools=expand(
            "results/QC/consensus/samtools/stats/{sample}.txt", sample=samples
        ),
        config_file=config["multiqc_config"],
    output:
        "results/QC/consensus/multiQC/multiqc_report.html",
    log:
        "logs/multiQC_consensus.log",
    conda:
        "../envs/multiqc.yaml"
    envmodules:
        "MultiQC/1.10.1",
    threads: cluster["multiqc"]["threads"]
    resources:
        mem_mb=cluster["multiqc"]["mem_mb"],
        runtime=cluster["multiqc"]["runtime"],
    shell:
        """
        multiqc \
            --config {input.config_file} \
            results/QC/consensus/qualimap/ \
            results/QC/consensus/samtools/ \
            -o results/QC/consensus/multiQC -f
        """


rule get_read_info:
    input:
        bam="results/mapped_bams/{sample}_mapped_merged.bam",
        script="workflow/scripts/get_read_info.py",
    output:
        "results/QC/read_info/{sample}_rinfo.txt",
    log:
        "logs/get_read_info_{sample}.log",
    conda:
        "../envs/get_read_info.yaml"
    threads: cluster["getreadinfo"]["threads"]
    resources:
        mem_mb=cluster["getreadinfo"]["mem_mb"],
        runtime=cluster["getreadinfo"]["runtime"],
    params:
        buffer_size=cluster["getreadinfo"]["buffer_size"],
    shell:
        """
        python {input.script} {input.bam} | \
            sort -k 2,6 -k 1,1 -u --parallel={threads} -S {params.buffer_size} -T {resources.tmpdir} | \
            cut -f 2-6 | \
            uniq --count > {output}
        """


rule format_read_info:
    input:
        "results/QC/read_info/{sample}_rinfo.txt",
    output:
        "results/QC/read_info/{sample}.txt.gz",
    log:
        "logs/format_read_info_{sample}.log",
    conda:
        "../envs/get_read_info.yaml"
    threads: cluster["formatreadinfo"]["threads"]
    resources:
        mem_mb=cluster["formatreadinfo"]["mem_mb"],
        runtime=cluster["formatreadinfo"]["runtime"],
    script:
        "../scripts/format_read_info.py"
