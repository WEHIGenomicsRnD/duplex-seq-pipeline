/*
 * PREPROCESSING subworkflow
 * FASTQ → unmapped BAM → extract/correct UMIs → align → merge BAM alignment
 */

include { PICARD_FASTQTOSAM                          } from '../../modules/local/picard/fastqtosam/main'
include { FGBIO_EXTRACT_UMIS                         } from '../../modules/local/fgbio/extract_umis/main'
include { FGBIO_CORRECT_UMIS                         } from '../../modules/local/fgbio/correct_umis/main'
include { PICARD_SAMTOFASTQ                          } from '../../modules/local/picard/samtofastq/main'
include { PICARD_MERGEBAMALIGNMENT                   } from '../../modules/local/picard/mergebamalignment/main'

workflow PREPROCESS {

    take:
    ch_reads       // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    ch_ref         // channel: path(ref)
    ch_dict        // channel: path(dict)
    ch_bwa_index   // channel: path(bwa_index)
    correct_umis   // val: boolean

    main:
    def ch_versions = Channel.empty()

    //
    // Convert FASTQ to unmapped BAM
    //
    PICARD_FASTQTOSAM(ch_reads)
    ch_versions = ch_versions.mix(PICARD_FASTQTOSAM.out.versions.first())

    //
    // Extract UMIs from unmapped BAM
    //
    FGBIO_EXTRACT_UMIS(PICARD_FASTQTOSAM.out.bam)
    ch_versions = ch_versions.mix(FGBIO_EXTRACT_UMIS.out.versions.first())

    //
    // Optionally correct UMIs
    //
    if (correct_umis) {
        def ch_umis = Channel.fromPath(params.umis, checkIfExists: true)
        FGBIO_CORRECT_UMIS(FGBIO_EXTRACT_UMIS.out.bam, ch_umis.collect())
        ch_versions = ch_versions.mix(FGBIO_CORRECT_UMIS.out.versions.first())
        ch_umi_bam = FGBIO_CORRECT_UMIS.out.bam
    } else {
        ch_umi_bam = FGBIO_EXTRACT_UMIS.out.bam
    }

    //
    // Convert BAM to interleaved FASTQ for alignment
    //
    PICARD_SAMTOFASTQ(ch_umi_bam)
    ch_versions = ch_versions.mix(PICARD_SAMTOFASTQ.out.versions.first())

    //
    // Align with BWA MEM (uses nf-core module)
    // Expects interleaved FASTQ with -p flag (set via ext.args)
    //
    def ch_align_input = PICARD_SAMTOFASTQ.out.fastq.map { meta, fastq ->
        [meta, fastq, true, false]   // [meta, reads, sort, index] – interleaved=true set via ext.args
    }

    // BWA_MEM is included from main workflow level
    // We emit the fastq for the parent workflow to wire through BWA_MEM

    //
    // Merge aligned BAM with unmapped BAM to preserve UMI tags
    //
    // NOTE: ch_mapped_bam is joined to ch_umi_bam in the parent workflow
    // after BWA_MEM runs, then passed to PICARD_MERGEBAMALIGNMENT

    emit:
    umi_bam       = ch_umi_bam                         // channel: [ val(meta), path(bam) ]
    fastq         = PICARD_SAMTOFASTQ.out.fastq        // channel: [ val(meta), path(fastq) ]
    versions      = ch_versions                         // channel: path(versions.yml)
}
