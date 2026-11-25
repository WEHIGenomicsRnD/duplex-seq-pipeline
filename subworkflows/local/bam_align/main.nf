//
// Realign an unmapped BAM: SamToFastq -> bwa mem -> MergeBamAlignment -> index.
// Used for the raw reads and, later, for the duplex consensus reads.
//

include { GATK4_SAMTOFASTQ } from '../../../modules/nf-core/gatk4/samtofastq/main'
include { BWA_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_MERGEBAMALIGNMENT } from '../../../modules/nf-core/gatk4/mergebamalignment/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_ALIGN {
    take:
    unmapped_bam // channel: [ val(meta), path(bam) ]
    fasta // channel: [ val(meta), path(fasta) ]
    dict // channel: [ val(meta), path(dict) ]
    bwa_index // channel: [ val(meta), path(index) ]

    main:

    GATK4_SAMTOFASTQ(unmapped_bam)

    // bwa mem reads the mates from two FASTQs; the unsorted BAM keeps FASTQ (query) order
    BWA_MEM(GATK4_SAMTOFASTQ.out.fastq, bwa_index, fasta, false)

    // MergeBamAlignment needs the aligned and unmapped BAMs in the same (query) order
    ch_merge_input = BWA_MEM.out.bam.join(unmapped_bam, failOnMismatch: true)

    GATK4_MERGEBAMALIGNMENT(ch_merge_input, fasta, dict)

    SAMTOOLS_INDEX(GATK4_MERGEBAMALIGNMENT.out.bam)

    ch_bam_bai = GATK4_MERGEBAMALIGNMENT.out.bam.join(SAMTOOLS_INDEX.out.index, failOnMismatch: true)

    emit:
    bam = GATK4_MERGEBAMALIGNMENT.out.bam // channel: [ val(meta), path(bam) ]
    bam_bai = ch_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
}
