//
// Variant calling on the filtered duplex consensus reads: clip overlapping read pairs, then call
// variants with bcftools mpileup and bin/filter_pileup.py (the default) or with the naive variant
// caller (nvc) and bin/filter_vcf.py. Both write the same VCF fields. With a target BED file only the
// positions in it are reported; nvc is then run on the reads that overlap it, as nvc is slow on
// off-target reads.
//

include { PICARD_SORTSAM as PICARD_SORTSAM_FILTERED } from '../../../modules/nf-core/picard/sortsam/main'
include { FGBIO_CLIPBAM } from '../../../modules/local/fgbio/clipbam/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_CLIPPED } from '../../../modules/nf-core/samtools/index/main'
include { BCFTOOLS_PILEUP } from '../../../modules/local/bcftools/pileup/main'
include { FILTER_PILEUP } from '../../../modules/local/filterpileup/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_TARGETS } from '../../../modules/nf-core/samtools/view/main'
include { NVC } from '../../../modules/local/nvc/main'
include { FILTER_VCF } from '../../../modules/local/filtervcf/main'

workflow VARIANT_CALLING {
    take:
    bam // channel: [ val(meta), path(bam) ] filtered consensus BAM, coordinate sorted
    fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]
    fasta_fai_dict // channel: [ val(meta), path(fasta), path(fai), path(dict) ]
    target_bed // value: path(bed) or [] (restricts the calls to the target regions)
    variant_caller // value: 'pileup' or 'nvc'

    main:

    // fgbio >= 2 needs query-grouped input to clip read pairs; ClipBam writes coordinate order
    PICARD_SORTSAM_FILTERED(bam, 'queryname')

    FGBIO_CLIPBAM(PICARD_SORTSAM_FILTERED.out.bam, fasta_fai_dict)

    SAMTOOLS_INDEX_CLIPPED(FGBIO_CLIPBAM.out.bam)

    ch_clipped_bam_bai = FGBIO_CLIPBAM.out.bam.join(SAMTOOLS_INDEX_CLIPPED.out.index, failOnMismatch: true)

    if (variant_caller == 'nvc') {
        // nvc's run time grows with the number of separately covered stretches of the genome, which
        // off-target reads (e.g. chimeric read pairs) multiply; with a target BED file it only gets
        // the reads that overlap the targets, and filter_vcf.py keeps the positions in them
        if (target_bed) {
            SAMTOOLS_VIEW_TARGETS(
                ch_clipped_bam_bai,
                [[:], [], []],
                [[:], []],
                [[:], target_bed],
                'bai',
            )
            ch_nvc_bam_bai = SAMTOOLS_VIEW_TARGETS.out.bam.join(SAMTOOLS_VIEW_TARGETS.out.bai, failOnMismatch: true)
        }
        else {
            ch_nvc_bam_bai = ch_clipped_bam_bai
        }

        NVC(ch_nvc_bam_bai, fasta_fai)

        FILTER_VCF(NVC.out.vcf, target_bed)

        ch_vcf = NVC.out.vcf
        ch_filtered_vcf = FILTER_VCF.out.vcf
    }
    else {
        BCFTOOLS_PILEUP(ch_clipped_bam_bai, fasta_fai, target_bed)

        FILTER_PILEUP(BCFTOOLS_PILEUP.out.vcf)

        ch_vcf = BCFTOOLS_PILEUP.out.vcf
        ch_filtered_vcf = FILTER_PILEUP.out.vcf
    }

    emit:
    clipped_bam_bai = ch_clipped_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    vcf = ch_vcf // channel: [ val(meta), path(vcf) ] nvc VCF or the mpileup records of the candidate positions
    filtered_vcf = ch_filtered_vcf // channel: [ val(meta), path(vcf) ]
}
