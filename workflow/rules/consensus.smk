rule group_reads_by_umi:
    input:
        "results/mapped_bams/{sample}_mapped_merged.bam",
    output:
        "results/grouped_by_umi/{sample}_grouped.bam",
    log:
        "logs/group_reads_by_umi_{sample}.log",
    conda:
        "../envs/fgbio.yaml"
    threads: cluster["fgbio"]["threads"]
    resources:
        mem_mb=cluster["fgbio"]["mem_mb"],
        runtime=cluster["fgbio"]["runtime"],
    params:
        min_mapq=config["min_mapq"],
    shell:
        """
        fgbio -Xmx{resources.mem_mb}m \
            -Djava.io.tmpdir={resources.tmpdir} \
            GroupReadsByUmi \
            --input={input} \
            --output={output} \
            --strategy=paired \
            --edits=1 \
            --min-map-q={params.min_mapq} \
        """


rule call_duplex_consensus:
    input:
        "results/grouped_by_umi/{sample}_grouped.bam",
    output:
        "results/consensus/{sample}_consensus_unmapped.bam",
    log:
        "logs/call_duplex_consensus_{sample}.log",
    conda:
        "../envs/fgbio.yaml"
    threads: cluster["fgbio"]["threads"]
    resources:
        mem_mb=cluster["fgbio"]["mem_mb"],
        runtime=cluster["fgbio"]["runtime"],
    params:
        error_rate_pre_umi=config["error_rate_pre_umi"],
        error_rate_post_umi=config["error_rate_post_umi"],
        min_input_baseq=config["min_input_baseq"],
    shell:
        """
        fgbio -Xmx{resources.mem_mb}m \
            -Djava.io.tmpdir={resources.tmpdir} \
            CallDuplexConsensusReads \
            --input={input} \
            --output={output} \
            --error-rate-pre-umi={params.error_rate_pre_umi} \
            --error-rate-post-umi={params.error_rate_post_umi} \
            --min-input-base-quality={params.min_input_baseq}
        """


rule remap_sam_to_fastq:
    input:
        "results/consensus/{sample}_consensus_unmapped.bam",
    output:
        "results/consensus/{sample}_consensus.fastq",
    log:
        "logs/remap_sam_to_fastq_{sample}.log",
    conda:
        "../envs/picard.yaml"
    threads: cluster["picard"]["threads"]
    resources:
        mem_mb=cluster["picard"]["mem_mb"],
        runtime=cluster["picard"]["runtime"],
    params:
        java_mem=cluster["picard"]["java_mem"],
    shell:
        """
        picard -Xmx{params.java_mem}m SamToFastq \
            I={input} \
            F={output} \
            INTERLEAVE=true \
            TMP_DIR={resources.tmpdir}
        """


rule remap_align:
    input:
        fastq="results/consensus/{sample}_consensus.fastq",
        ref=config["ref"],
    output:
        "results/consensus/{sample}_consensus_mapped.bam",
    log:
        "logs/remap_align_{sample}.log",
    conda:
        "../envs/bwa_samtools.yaml"
    envmodules:
        "bwa/0.7.17",
    threads: cluster["bwa"]["threads"]
    resources:
        mem_mb=cluster["bwa"]["mem_mb"],
        runtime=cluster["bwa"]["runtime"],
    shell:
        """
        bwa mem -p -t {threads} {input.ref} {input.fastq} | samtools view -bS - > {output}
        """


rule remap_merge_bam_alignment:
    input:
        unmapped="results/consensus/{sample}_consensus_unmapped.bam",
        mapped="results/consensus/{sample}_consensus_mapped.bam",
        ref=config["ref"],
    output:
        "results/consensus/{sample}_consensus_mapped_merged.bam",
    log:
        "logs/remap_merge_bam_alignment_{sample}.log",
    conda:
        "../envs/picard.yaml"
    threads: cluster["picard"]["threads"]
    resources:
        mem_mb=cluster["picard"]["mem_mb"],
        runtime=cluster["picard"]["runtime"],
    params:
        java_mem=cluster["picard"]["java_mem"],
    shell:
        """
        picard -Xmx{params.java_mem}m MergeBamAlignment \
                UNMAPPED={input.unmapped} \
                ALIGNED={input.mapped} \
                O={output} \
                R={input.ref} \
                SO=coordinate \
                ALIGNER_PROPER_PAIR_FLAGS=true \
                MAX_GAPS=-1 \
                ORIENTATIONS=FR \
                CREATE_INDEX=true \
                TMP_DIR={resources.tmpdir}
        """


rule filter_consensus_reads:
    input:
        bam="results/consensus/{sample}_consensus_mapped_merged.bam",
        ref=config["ref"],
    output:
        "results/consensus/{sample}_consensus_mapped_merged_filtered.bam",
    log:
        "logs/filter_consensus_reads_{sample}.log",
    conda:
        "../envs/fgbio.yaml"
    threads: cluster["fgbio"]["threads"]
    resources:
        mem_mb=cluster["fgbio"]["mem_mb"],
        runtime=cluster["fgbio"]["runtime"],
    params:
        min_reads=config["min_reads"],
        max_read_error_rate=config["max_read_error_rate"],
        max_base_error_rate=config["max_base_error_rate"],
        min_base_quality=config["min_base_quality"],
        max_no_call_fraction=config["max_no_call_fraction"],
    shell:
        """
        fgbio -Xmx{resources.mem_mb}m \
            -Djava.io.tmpdir={resources.tmpdir} \
            FilterConsensusReads \
            --input={input.bam} \
            --output={output} \
            --ref {input.ref} \
            --min-reads={params.min_reads} \
            --max-read-error-rate={params.max_read_error_rate} \
            --max-base-error-rate={params.max_base_error_rate} \
            --min-base-quality={params.min_base_quality} \
            --max-no-call-fraction={params.max_no_call_fraction}
        """
