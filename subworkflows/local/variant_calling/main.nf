//
// Variant calling on the filtered duplex consensus reads: clip overlapping read pairs, call
// variants with the naive variant caller (nvc) and filter the VCF.
//

include { PICARD_SORTSAM as PICARD_SORTSAM_FILTERED } from '../../../modules/nf-core/picard/sortsam/main'
include { FGBIO_CLIPBAM } from '../../../modules/local/fgbio/clipbam/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_CLIPPED } from '../../../modules/nf-core/samtools/index/main'
include { NVC } from '../../../modules/local/nvc/main'
include { FILTER_VCF } from '../../../modules/local/filtervcf/main'

workflow VARIANT_CALLING {
    take:
    bam // channel: [ val(meta), path(bam) ] filtered consensus BAM, coordinate sorted
    fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]
    fasta_fai_dict // channel: [ val(meta), path(fasta), path(fai), path(dict) ]

    main:

    // fgbio >= 2 needs query-grouped input to clip read pairs; ClipBam writes coordinate order
    PICARD_SORTSAM_FILTERED(bam, 'queryname')

    FGBIO_CLIPBAM(PICARD_SORTSAM_FILTERED.out.bam, fasta_fai_dict)

    SAMTOOLS_INDEX_CLIPPED(FGBIO_CLIPBAM.out.bam)

    ch_clipped_bam_bai = FGBIO_CLIPBAM.out.bam.join(SAMTOOLS_INDEX_CLIPPED.out.index, failOnMismatch: true)

    NVC(ch_clipped_bam_bai, fasta_fai)

    FILTER_VCF(NVC.out.vcf)

    emit:
    clipped_bam_bai = ch_clipped_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    vcf = NVC.out.vcf // channel: [ val(meta), path(vcf) ]
    filtered_vcf = FILTER_VCF.out.vcf // channel: [ val(meta), path(vcf) ]
}
