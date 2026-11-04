/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/gunzip/main'
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference/main'
include { FASTQ_ALIGN_UMI } from '../subworkflows/local/fastq_align_umi/main'
include { CONSENSUS } from '../subworkflows/local/consensus/main'
include { QC } from '../subworkflows/local/qc/main'
include { VARIANT_CALLING } from '../subworkflows/local/variant_calling/main'
include { softwareVersionsToYAML } from '../subworkflows/local/utils_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DUPLEXSEQ {
    take:
    ch_samplesheet // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ], path(umi_file) or [] ]

    main:

    // reference files; indexes that are not supplied are built
    def ref_meta = [id: file(params.fasta.replaceFirst(/\.gz$/, '')).baseName]
    ch_fasta_input = channel.fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> [ref_meta, fasta] }
        .branch { _meta, fasta ->
            gzipped: fasta.name.endsWith('.gz')
            plain: true
        }

    // a gzipped reference is decompressed (samtools faidx and bwa need plain or bgzip FASTA)
    GUNZIP_FASTA(ch_fasta_input.gzipped)
    ch_fasta = ch_fasta_input.plain.mix(GUNZIP_FASTA.out.gunzip)
    ch_fasta_fai = params.fasta_fai ? channel.fromPath(params.fasta_fai, checkIfExists: true).map { fai -> [ref_meta, fai] } : channel.empty()
    ch_dict = params.dict ? channel.fromPath(params.dict, checkIfExists: true).map { dict -> [ref_meta, dict] } : channel.empty()
    ch_bwa_index = params.bwa ? channel.fromPath(params.bwa, checkIfExists: true, type: 'dir').map { index -> [ref_meta, index] } : channel.empty()

    PREPARE_REFERENCE(ch_fasta, ch_fasta_fai, ch_dict, ch_bwa_index)

    // FASTQ -> UMI-extracted, aligned BAM
    ch_reads = ch_samplesheet.map { meta, reads, _umis -> [meta, reads] }
    ch_umis = ch_samplesheet.map { meta, _reads, umis -> [meta, umis] }

    ch_fasta_value = ch_fasta.first()

    FASTQ_ALIGN_UMI(
        ch_reads,
        ch_umis,
        ch_fasta_value,
        PREPARE_REFERENCE.out.dict,
        PREPARE_REFERENCE.out.bwa_index,
    )

    // fasta bundled with its indexes, as expected by the fgbio, picard and samtools modules
    ch_fasta_fai = ch_fasta_value
        .combine(PREPARE_REFERENCE.out.fai)
        .map { meta, fasta, _meta_fai, fai -> [meta, fasta, fai] }
        .first()
    ch_fasta_fai_dict = ch_fasta_value
        .combine(PREPARE_REFERENCE.out.fai)
        .combine(PREPARE_REFERENCE.out.dict)
        .map { meta, fasta, _meta_fai, fai, _meta_dict, dict -> [meta, fasta, fai, dict] }
        .first()

    // duplex consensus calling, realignment and filtering
    CONSENSUS(
        FASTQ_ALIGN_UMI.out.bam,
        ch_fasta_value,
        ch_fasta_fai_dict,
        PREPARE_REFERENCE.out.dict,
        PREPARE_REFERENCE.out.bwa_index,
        params.min_input_baseq,
        params.min_reads,
        params.min_base_quality,
        params.max_base_error_rate,
    )

    // variant calling on the filtered consensus reads
    VARIANT_CALLING(
        CONSENSUS.out.filtered_bam_bai.map { meta, bam, _bai -> [meta, bam] },
        ch_fasta_fai,
        ch_fasta_fai_dict,
    )

    // QC and a single MultiQC report; FastQ Screen only runs when a database folder is given
    if (params.fastq_screen_db && !file("${params.fastq_screen_db}/fastq_screen.conf").exists()) {
        error("--fastq_screen_db must be a folder containing 'fastq_screen.conf': ${params.fastq_screen_db}")
    }
    ch_fastq_screen_db = params.fastq_screen_db ? channel.value(file(params.fastq_screen_db, checkIfExists: true, type: 'dir')) : channel.empty()
    // the packaged MultiQC config is always used; --multiqc_config is added on top of it
    ch_multiqc_config = params.multiqc_config
        ? [file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true), file(params.multiqc_config, checkIfExists: true)]
        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    ch_target_bed = params.target_bed ? file(params.target_bed, checkIfExists: true) : []

    QC(
        ch_reads,
        FASTQ_ALIGN_UMI.out.bam_bai,
        CONSENSUS.out.filtered_bam_bai,
        ch_fasta_fai,
        ch_target_bed,
        ch_fastq_screen_db,
        params.run_qualimap,
        params.output_read_info,
        ch_multiqc_config,
    )

    // software versions from the `versions` topic
    ch_versions = softwareVersionsToYAML(channel.topic('versions'))

    emit:
    bam = FASTQ_ALIGN_UMI.out.bam // channel: [ val(meta), path(bam) ]
    consensus_bam = CONSENSUS.out.filtered_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    vcf = VARIANT_CALLING.out.filtered_vcf // channel: [ val(meta), path(vcf) ]
    multiqc_report = QC.out.report // channel: [ val(meta), path(html) ]
    versions = ch_versions // channel: path(software_versions.yml)
}
