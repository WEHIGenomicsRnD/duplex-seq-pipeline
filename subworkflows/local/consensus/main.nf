//
// Duplex consensus: group by UMI -> call duplex consensus -> realign -> filter -> index.
//

include { FGBIO_GROUPREADSBYUMI } from '../../../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_CALLDUPLEXCONSENSUSREADS } from '../../../modules/nf-core/fgbio/callduplexconsensusreads/main'
include { PICARD_SORTSAM } from '../../../modules/nf-core/picard/sortsam/main'
include { PICARD_SORTSAM as PICARD_SORTSAM_MAPPED } from '../../../modules/nf-core/picard/sortsam/main'
include { FGBIO_FILTERCONSENSUSREADS } from '../../../modules/nf-core/fgbio/filterconsensusreads/main'
include { BAM_ALIGN as BAM_ALIGN_CONSENSUS } from '../bam_align/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED } from '../../../modules/nf-core/samtools/index/main'

workflow CONSENSUS {
    take:
    bam // channel: [ val(meta), path(bam) ] aligned, UMI-tagged BAM
    fasta // channel: [ val(meta), path(fasta) ]
    fasta_fai_dict // channel: [ val(meta), path(fasta), path(fai), path(dict) ]
    dict // channel: [ val(meta), path(dict) ]
    bwa_index // channel: [ val(meta), path(index) ]
    min_input_baseq // value: minimum input base quality for consensus calling
    min_reads // value: FilterConsensusReads minimum reads ("duplex AB BA")
    min_base_quality // value: FilterConsensusReads minimum base quality
    max_base_error_rate // value: FilterConsensusReads maximum base error rate

    main:

    // strategy "paired" is what duplex calling needs
    FGBIO_GROUPREADSBYUMI(bam, 'paired')

    // Snakemake keeps the fgbio default of one read per strand here and filters afterwards
    FGBIO_CALLDUPLEXCONSENSUSREADS(FGBIO_GROUPREADSBYUMI.out.bam, 1, min_input_baseq)

    // fgbio names consensus reads numerically (:1, :2, ... :10) but MergeBamAlignment needs the
    // unmapped BAM in Picard's lexicographic queryname order, to match the aligned reads
    PICARD_SORTSAM(FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam, 'queryname')

    BAM_ALIGN_CONSENSUS(PICARD_SORTSAM.out.bam, fasta, dict, bwa_index)

    // fgbio >= 2 requires query-grouped input for FilterConsensusReads; it writes coordinate order
    PICARD_SORTSAM_MAPPED(BAM_ALIGN_CONSENSUS.out.bam, 'queryname')

    FGBIO_FILTERCONSENSUSREADS(
        PICARD_SORTSAM_MAPPED.out.bam,
        fasta_fai_dict,
        min_reads,
        min_base_quality,
        max_base_error_rate,
    )

    SAMTOOLS_INDEX_FILTERED(FGBIO_FILTERCONSENSUSREADS.out.bam)

    ch_filtered_bam_bai = FGBIO_FILTERCONSENSUSREADS.out.bam.join(SAMTOOLS_INDEX_FILTERED.out.index, failOnMismatch: true)

    emit:
    unmapped_bam = FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam // channel: [ val(meta), path(bam) ]
    mapped_bam = BAM_ALIGN_CONSENSUS.out.bam // channel: [ val(meta), path(bam) ]
    filtered_bam_bai = ch_filtered_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    grouping_metrics = FGBIO_GROUPREADSBYUMI.out.read_metrics // channel: [ val(meta), path(txt) ]
    family_histogram = FGBIO_GROUPREADSBYUMI.out.histogram // channel: [ val(meta), path(txt) ]
}
