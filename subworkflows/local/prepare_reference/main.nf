//
// Prepare the reference: samtools faidx, picard CreateSequenceDictionary and bwa index.
// Each index is only built when it is not supplied.
//

include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { BWA_INDEX } from '../../../modules/nf-core/bwa/index/main'

workflow PREPARE_REFERENCE {
    take:
    fasta // channel: [ val(meta), path(fasta) ]
    fasta_fai // channel: [ val(meta), path(fai) ] or empty
    dict // channel: [ val(meta), path(dict) ] or empty
    bwa_index // channel: [ val(meta), path(index) ] or empty

    main:

    // samtools faidx
    ch_fasta_for_faidx = fasta
        .join(fasta_fai, remainder: true)
        .map { meta, fa, fai -> [meta, fa, fai ?: []] }
        .branch { _meta, _fa, fai ->
            have: fai
            missing: true
        }

    SAMTOOLS_FAIDX(ch_fasta_for_faidx.missing, false)
    ch_fai = ch_fasta_for_faidx.have
        .map { meta, _fa, fai -> [meta, fai] }
        .mix(SAMTOOLS_FAIDX.out.fai)

    // picard CreateSequenceDictionary
    ch_fasta_for_dict = fasta
        .join(dict, remainder: true)
        .branch { _meta, _fa, d ->
            have: d
            missing: true
        }

    PICARD_CREATESEQUENCEDICTIONARY(ch_fasta_for_dict.missing.map { meta, fa, _d -> [meta, fa] })
    ch_dict = ch_fasta_for_dict.have
        .map { meta, _fa, d -> [meta, d] }
        .mix(PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict)

    // bwa index
    ch_fasta_for_bwa = fasta
        .join(bwa_index, remainder: true)
        .branch { _meta, _fa, idx ->
            have: idx
            missing: true
        }

    BWA_INDEX(ch_fasta_for_bwa.missing.map { meta, fa, _idx -> [meta, fa] })
    ch_bwa_index = ch_fasta_for_bwa.have
        .map { meta, _fa, idx -> [meta, idx] }
        .mix(BWA_INDEX.out.index)

    emit:
    fai = ch_fai.first() // value channel: [ val(meta), path(fai) ]
    dict = ch_dict.first() // value channel: [ val(meta), path(dict) ]
    bwa_index = ch_bwa_index.first() // value channel: [ val(meta), path(index) ]
}
