/*
 * CONSENSUS subworkflow
 * Group reads by UMI → call duplex consensus → remap → filter
 */

include { PICARD_SAMTOFASTQ as PICARD_SAMTOFASTQ_CONSENSUS        } from '../../modules/local/picard/samtofastq/main'
include { PICARD_MERGEBAMALIGNMENT as PICARD_MERGEBAMALIGNMENT_CONSENSUS } from '../../modules/local/picard/mergebamalignment/main'

workflow CONSENSUS {

    take:
    ch_merged_bam  // channel: [ val(meta), path(bam) ] - merged mapped BAM
    ch_ref         // channel: path(ref)
    ch_dict        // channel: path(dict)
    ch_bwa_index   // channel: path(bwa_index)

    main:
    def ch_versions = Channel.empty()

    //
    // Group reads by UMI (nf-core module)
    // Wired from parent workflow: FGBIO_GROUPREADSBYUMI
    //

    //
    // Call duplex consensus reads (nf-core module)
    // Wired from parent workflow: FGBIO_CALLDUPLEXCONSENSUSREADS
    //

    //
    // Convert consensus unmapped BAM to interleaved FASTQ for remapping
    //
    PICARD_SAMTOFASTQ_CONSENSUS(ch_merged_bam)
    ch_versions = ch_versions.mix(PICARD_SAMTOFASTQ_CONSENSUS.out.versions.first())

    emit:
    fastq    = PICARD_SAMTOFASTQ_CONSENSUS.out.fastq   // channel: [ val(meta), path(fastq) ]
    versions = ch_versions                               // channel: path(versions.yml)
}
