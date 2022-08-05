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
    threads: cluster["fgbio"]["threads"]
    resources:
        mem_mb=cluster["fgbio"]["mem_mb"],
        runtime=cluster["fgbio"]["runtime"],
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
    threads: cluster["varscan"]["threads"]
    resources:
        mem_mb=cluster["varscan"]["mem_mb"],
        runtime=cluster["varscan"]["runtime"],
    params:
        min_coverage=config["min_coverage"],
        min_reads2=config["min_reads2"],
        min_vaf=config["min_vaf"],
    shell:
        """
        samtools mpileup -f {input.ref} {input.bam} | \
            varscan mpileup2snp -Xmx{resources.mem_mb} \
            --min-avg-qual 0 \
            --min-coverage {params.min_coverage} \
            --min-reads2 {params.min_reads2} \
            --min-var-freq {params.min_vaf} \
            --strand-filter 0 \
            --p-value 1 \
            --output-vcf 1 > {output}
        """
