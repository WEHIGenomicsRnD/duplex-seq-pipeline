/*
 * VARIANT_CALLING subworkflow
 * Clip BAM → call variants (VarScan or NVC) → filter VCF
 */

include { FGBIO_CLIPBAM } from '../../modules/local/fgbio/clipbam/main'
include { VARSCAN       } from '../../modules/local/varscan/main'
include { NVC           } from '../../modules/local/nvc/main'
include { FILTER_VCF    } from '../../modules/local/filter_vcf/main'

workflow VARIANT_CALLING {

    take:
    ch_filtered_bam   // channel: [ val(meta), path(bam) ]
    ch_ref            // channel: path(ref)
    ch_dict           // channel: path(dict)
    variant_caller    // val: 'varscan' or 'nvc'

    main:
    def ch_versions = Channel.empty()

    //
    // Clip overlapping reads
    //
    FGBIO_CLIPBAM(ch_filtered_bam, ch_ref, ch_dict)
    ch_versions = ch_versions.mix(FGBIO_CLIPBAM.out.versions.first())

    //
    // Call variants
    //
    if (variant_caller == 'varscan') {
        VARSCAN(FGBIO_CLIPBAM.out.bam, ch_ref)
        ch_versions = ch_versions.mix(VARSCAN.out.versions.first())
        ch_vcf = VARSCAN.out.vcf
    } else {
        NVC(FGBIO_CLIPBAM.out.bam, ch_ref)
        ch_versions = ch_versions.mix(NVC.out.versions.first())
        ch_vcf = NVC.out.vcf
    }

    //
    // Filter VCF (add DP field, remove N alleles)
    //
    FILTER_VCF(ch_vcf)
    ch_versions = ch_versions.mix(FILTER_VCF.out.versions.first())

    emit:
    vcf      = FILTER_VCF.out.vcf    // channel: [ val(meta), path(vcf) ]
    versions = ch_versions            // channel: path(versions.yml)
}
