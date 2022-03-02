rule create_dict:
    input:
        config["ref"],
    output:
        get_create_dict_output(),
    log:
        "logs/create_dict.log",
    conda:
        "../envs/picard.yaml"
    threads: 1
    resources:
        mem_mb=32768,
        runtime="1-0:0:0",
    shell:
        """
        picard -Xmx{resources.mem_mb}m CreateSequenceDictionary \
            R={input} \
            O={output} \
            TMP_DIR={resources.tmpdir}
        """


rule convert_to_unmapped_bam:
    input:
        r1="fastq/{sample}_R1.fastq.gz",
        r2="fastq/{sample}_R2.fastq.gz",
    output:
        "results/unmapped_bams/{sample}_unmapped.bam",
    log:
        "logs/convert_to_unmapped_bam_{sample}.log",
    conda:
        "../envs/picard.yaml"
    threads: 1
    resources:
        mem_mb=32768,
        runtime="1-0:0:0",
    shell:
        """
        picard -Xmx{resources.mem_mb}m FastqToSam \
            F1={input.r1} \
            F2={input.r2} \
            O={output} \
            SM={wildcards.sample} \
            RG=null \
            TMP_DIR={resources.tmpdir}
        """


rule extract_umis:
    input:
        "results/unmapped_bams/{sample}_unmapped.bam",
    output:
        "results/unmapped_bams/{sample}_unmapped_{umi_postfix}wUMI.bam".format(
            sample="{sample}", umi_postfix=umi_postfix
        ),
    log:
        "logs/extract_umis_{sample}.log",
    conda:
        "../envs/fgbio.yaml"
    threads: 1
    resources:
        mem_mb=32768,
        runtime="0-12:0:0",
    params:
        read_structure=config["read_structure"],
        molecular_index_tags=config["molecular_index_tags"],
        single_tag=config["single_tag"],
    shell:
        """
        fgbio -Xmx{resources.mem_mb}m \
            -Djava.io.tmpdir={resources.tmpdir} \
            ExtractUmisFromBam \
            --input={input} \
            --output={output} \
            --read-structure={params.read_structure} \
            --molecular-index-tags={params.molecular_index_tags} \
            --single-tag={params.single_tag}
        """


if correct_umis:

    rule correct_umis:
        input:
            bam="results/unmapped_bams/{sample}_unmapped_uncorrected_wUMI.bam",
            umis=config["umis"],
        output:
            "results/unmapped_bams/{sample}_unmapped_wUMI.bam",
        log:
            "logs/correct_umis_{sample}.log",
        conda:
            "../envs/fgbio.yaml"
        threads: 1
        resources:
            mem_mb=32768,
            runtime="0-12:0:0",
        params:
            max_mismatches=config["max_mismatches"],
            min_distance=config["min_distance"],
        shell:
            """
            fgbio -Xmx{resources.mem_mb}m \
                -Djava.io.tmpdir={resources.tmpdir} \
                CorrectUmis \
                --input={input.bam} \
                --output={output} \
                --max-mismatches={params.max_mismatches} \
                --min-distance={params.min_distance} \
                -t RX -U {input.umis}
            """


rule bam_to_fastq:
    input:
        "results/unmapped_bams/{sample}_unmapped_wUMI.bam",
    output:
        "results/unmapped_bams/{sample}.fastq",
    log:
        "logs/sam_to_fastq_{sample}.log",
    conda:
        "../envs/picard.yaml"
    threads: 1
    resources:
        mem_mb=32768,
        runtime="0-12:0:0",
    shell:
        """
        picard -Xmx{resources.mem_mb}m SamToFastq \
            I={input} \
            F={output} \
            INTERLEAVE=true \
            TMP_DIR={resources.tmpdir}
        """


rule align:
    input:
        fastq="results/unmapped_bams/{sample}.fastq",
        ref=config["ref"],
    output:
        "results/mapped_bams/{sample}_mapped.bam",
    log:
        "logs/align_{sample}.log",
    conda:
        "../envs/bwa_samtools.yaml"
    envmodules:
        "bwa/0.7.17",
    threads: 16
    resources:
        mem_mb=36864,
        runtime="0-12:0:0",
    shell:
        """
        bwa mem -p -t {threads} {input.ref} {input.fastq} | samtools view -bS - > {output}
        """


rule merge_bam_alignment:
    input:
        unmapped="results/unmapped_bams/{sample}_unmapped_wUMI.bam",
        mapped="results/mapped_bams/{sample}_mapped.bam",
        ref=config["ref"],
    output:
        "results/mapped_bams/{sample}_mapped_merged.bam",
    log:
        "logs/merge_bam_alignment_{sample}.log",
    conda:
        "../envs/picard.yaml"
    threads: 1
    resources:
        mem_mb=32768,
        runtime="0-12:0:0",
    shell:
        """
        picard -Xmx{resources.mem_mb}m MergeBamAlignment \
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
    threads: 1
    resources:
        mem_mb=32768,
        runtime="0-12:0:0",
    shell:
        """
        picard -Xmx{resources.mem_mb}m MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            TMP_DIR={resources.tmpdir}
        """
