//
// Merge the lanes of the samples that have several (see groupLanes in utils_pipeline): the files of the
// lanes are concatenated, in lane order, into one pair of FASTQ files per sample. Samples with one lane
// pass through untouched.
//

include { CAT_FASTQ } from '../../../modules/nf-core/cat/fastq/main'

workflow MERGE_LANES {
    take:
    samplesheet // channel: [ val(meta), [ path(fastq) ... ], path(umi_file) or [] ]; meta.n_lanes is the number of lanes

    main:

    ch_by_lanes = samplesheet.branch { meta, _reads, _umis ->
        several: meta.n_lanes > 1
        single: true
    }

    // CAT_FASTQ expects the files of the lanes one after the other: R1, R2, R1, R2, ...
    CAT_FASTQ(ch_by_lanes.several.map { meta, reads, _umis -> [meta, reads] })

    ch_merged = CAT_FASTQ.out.reads.join(ch_by_lanes.several.map { meta, _reads, umis -> [meta, umis] }, failOnMismatch: true)

    // the number of lanes is only needed here: afterwards meta is the same as for a run without merging
    ch_samples = ch_by_lanes.single
        .mix(ch_merged)
        .map { meta, reads, umis -> [meta.findAll { key, _value -> key != 'n_lanes' }, reads, umis] }

    emit:
    samplesheet = ch_samples // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ], path(umi_file) or [] ]
}
