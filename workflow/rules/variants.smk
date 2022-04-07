rule clip_bam:
    input:
        bam="results/consensus/{sample}_consensus_mapped_merged_filtered.bam",
        ref=config["ref"],
    output:
        "results/consensus/{sample}_consensus_mapped_merged_filtered_clipped.bam",
    log:
        "logs/clip_bam_{sample}.log",
    conda:
        "../envs/fgbio.yaml"
    threads: 1
    resources:
        mem_mb=32768,
        runtime="0-12:0:0",
    shell:
        """
        fgbio -Xmx{resources.mem_mb}m \
            -Djava.io.tmpdir={resources.tmpdir} \
            ClipBam \
            --input={input.bam} \
            --output={output} \
            --ref {input.ref} \
            --clipping-mode=Hard \
            --clip-overlapping-reads=true
        """


rule call_variants:
    input:
        bam="results/consensus/{sample}_consensus_mapped_merged_filtered_clipped.bam",
        ref=config["ref"],
    output:
        "results/variants/{sample}.vcf",
    log:
        "logs/varscan_{sample}.log",
    envmodules:
        "samtools/1.9",
        "varscan/2.3.9",
    conda:
        "../envs/varscan.yaml"
    threads: 1
    resources:
        mem_mb=32768,
        runtime="1-0:0:0",
    params:
        min_coverage=config["min_coverage"],
        min_reads=config["min_reads"],
        min_vaf=config["min_vaf"],
        strand_filter=config["strand_filter"],
    shell:
        """
        samtools mpileup -f {input.ref} {input.bam} | \
            varscan mpileup2snp -Xmx{resources.mem_mb} \
            --min-avg-qual 0 \
            --min-coverage {params.min_coverage} \
            --min-reads2 {params.min_reads} \
            --min-var-freq {params.min_vaf} \
            --strand-filter {params.strand_filter} \
            --p-value 1 \
            --output-vcf 1 > {output}
        """
