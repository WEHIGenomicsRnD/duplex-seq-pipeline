/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/gunzip/main'
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference/main'
include { FASTQ_ALIGN_UMI } from '../subworkflows/local/fastq_align_umi/main'
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

    FASTQ_ALIGN_UMI(
        ch_reads,
        ch_umis,
        ch_fasta.first(),
        PREPARE_REFERENCE.out.dict,
        PREPARE_REFERENCE.out.bwa_index,
    )

    // software versions from the `versions` topic
    ch_versions = softwareVersionsToYAML(channel.topic('versions'))

    emit:
    bam = FASTQ_ALIGN_UMI.out.bam // channel: [ val(meta), path(bam) ]
    versions = ch_versions // channel: path(software_versions.yml)
}
