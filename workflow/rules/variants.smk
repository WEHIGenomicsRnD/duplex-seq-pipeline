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


if variant_caller == "varscan":

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
                varscan mpileup2snp -Xmx{resources.mem_mb}m \
                --min-avg-qual 0 \
                --min-coverage {params.min_coverage} \
                --min-reads2 {params.min_reads2} \
                --min-var-freq {params.min_vaf} \
                --strand-filter 0 \
                --p-value 1 \
                --output-vcf 1 > {output}
            """


elif variant_caller == "nvc":

    rule call_variants:
        input:
            bam="results/consensus/{sample}_consensus_mapped_merged_filtered_clipped.bam",
            ref=config["ref"],
        output:
            "results/variants/{sample}.vcf",
        log:
            "logs/nvc_{sample}.log",
        conda:
            "../envs/nvc.yaml"
        threads: 1
        resources:
            mem_mb=cluster["varscan"]["mem_mb"],
            runtime=cluster["varscan"]["runtime"],
        params:
            ploidy=config["ploidy"],
            min_reads2=config["min_reads2"],
            min_mapq=config["min_mapq"],
        shell:
            """
            samtools index {input.bam} ;
            naive_variant_caller.py \
                --bam={input.bam} \
                --output_vcf_filename={output} \
                --reference_genome_filename={input.ref} \
                --variants_only \
                --ploidy={params.ploidy} \
                --min_support_depth={params.min_reads2} \
                --min_mapping_quality={params.min_mapq} \
                --min_base_quality=0
            """


    rule filter_vcf:
       input:
           "results/variants/{sample}.vcf",
       output:
           "results/variants/{sample}_filtered.vcf"",
       log:
           "logs/filter_vcf_{sample}.log",
       conda:
           "../envs/get_read_info.yaml"
       threads: 1
       resources:
           mem_mb=cluster["varscan"]["mem_mb"],
           runtime=cluster["varscan"]["runtime"],
       script:
           "workflow/scripts/filter_vcf.py"
