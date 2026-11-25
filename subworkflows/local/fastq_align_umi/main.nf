//
// FASTQ -> unmapped BAM -> extract (and optionally correct) UMIs -> aligned BAM.
//

include { PICARD_FASTQTOSAM } from '../../../modules/nf-core/picard/fastqtosam/main'
include { FGBIO_EXTRACTUMISFROMBAM } from '../../../modules/local/fgbio/extractumisfrombam/main'
include { FGBIO_CORRECTUMIS } from '../../../modules/local/fgbio/correctumis/main'
include { BAM_ALIGN } from '../bam_align/main'

workflow FASTQ_ALIGN_UMI {
    take:
    reads // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ] ]
    umis // channel: [ val(meta), path(umi_file) or [] ]
    fasta // channel: [ val(meta), path(fasta) ]
    dict // channel: [ val(meta), path(dict) ]
    bwa_index // channel: [ val(meta), path(index) ]

    main:

    PICARD_FASTQTOSAM(reads)

    FGBIO_EXTRACTUMISFROMBAM(PICARD_FASTQTOSAM.out.bam)

    // only correct UMIs for samples that have a UMI file
    ch_extracted = FGBIO_EXTRACTUMISFROMBAM.out.bam
        .join(umis, failOnMismatch: true)
        .branch { _meta, _bam, umi_file ->
            correct: umi_file
            skip: true
        }

    FGBIO_CORRECTUMIS(ch_extracted.correct)

    ch_unmapped_umi = ch_extracted.skip
        .map { meta, bam, _umi_file -> [meta, bam] }
        .mix(FGBIO_CORRECTUMIS.out.bam)

    BAM_ALIGN(ch_unmapped_umi, fasta, dict, bwa_index)

    emit:
    unmapped_bam = ch_unmapped_umi // channel: [ val(meta), path(bam) ]
    bam = BAM_ALIGN.out.bam // channel: [ val(meta), path(bam) ]
    bam_bai = BAM_ALIGN.out.bam_bai // channel: [ val(meta), path(bam), path(bai) ]
}
