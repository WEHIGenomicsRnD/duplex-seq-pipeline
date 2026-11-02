//
// QC: FastQC, FastQ Screen, MarkDuplicates metrics, samtools stats/flagstat, Qualimap and the duplex
// metrics on the mapped and consensus BAMs, summarised in a single MultiQC report.
// Consensus outputs get a `_consensus` prefix so they are separate samples in MultiQC.
//

include { FASTQC } from '../../../modules/nf-core/fastqc/main'
include { FASTQSCREEN_FASTQSCREEN as FASTQSCREEN } from '../../../modules/nf-core/fastqscreen/fastqscreen/main'
include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { QUALIMAP_BAMQC } from '../../../modules/nf-core/qualimap/bamqc/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_CONSENSUS } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_CONSENSUS } from '../../../modules/nf-core/samtools/flagstat/main'
include { QUALIMAP_BAMQC as QUALIMAP_BAMQC_CONSENSUS } from '../../../modules/nf-core/qualimap/bamqc/main'
include { CALC_DUPLEX_METRICS } from '../../../modules/local/calcduplexmetrics/main'
include { MULTIQC } from '../../../modules/nf-core/multiqc/main'

workflow QC {
    take:
    reads // channel: [ val(meta), [ path(fastq_1), path(fastq_2) ] ]
    bam_bai // channel: [ val(meta), path(bam), path(bai) ] mapped, UMI-tagged BAM
    consensus_bam_bai // channel: [ val(meta), path(bam), path(bai) ] filtered consensus BAM
    fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]
    target_bed // value: path(bed) or [] (enables the on-target duplex metrics)
    fastq_screen_db // channel: path(dir with fastq_screen.conf), or empty to skip FastQ Screen
    run_qualimap // value: boolean
    multiqc_config // value: path(config) or [path(config), path(config)...]

    main:

    ch_multiqc_files = channel.empty()

    // raw reads
    FASTQC(reads)
    FASTQSCREEN(reads, fastq_screen_db)

    // mapped BAM
    PICARD_MARKDUPLICATES(bam_bai.map { meta, bam, _bai -> [meta, bam] }, fasta_fai)
    SAMTOOLS_STATS(bam_bai, fasta_fai)
    SAMTOOLS_FLAGSTAT(bam_bai)

    // duplex metrics are computed from the mapped BAM (as in Snakemake), not the consensus BAM
    CALC_DUPLEX_METRICS(
        bam_bai.map { meta, bam, _bai -> [meta, bam] },
        fasta_fai.map { meta, fasta, _fai -> [meta, fasta] },
        target_bed,
    )

    // consensus BAM
    SAMTOOLS_STATS_CONSENSUS(consensus_bam_bai, fasta_fai)
    SAMTOOLS_FLAGSTAT_CONSENSUS(consensus_bam_bai)

    ch_multiqc_files = ch_multiqc_files
        .mix(FASTQC.out.zip.map { _meta, f -> f })
        .mix(FASTQSCREEN.out.txt.map { _meta, f -> f })
        .mix(PICARD_MARKDUPLICATES.out.metrics.map { _meta, f -> f })
        .mix(SAMTOOLS_STATS.out.stats.map { _meta, f -> f })
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { _meta, f -> f })
        .mix(CALC_DUPLEX_METRICS.out.csv.map { _meta, f -> f })
        .mix(SAMTOOLS_STATS_CONSENSUS.out.stats.map { _meta, f -> f })
        .mix(SAMTOOLS_FLAGSTAT_CONSENSUS.out.flagstat.map { _meta, f -> f })

    // Qualimap is resource-intensive, so it is opt-in (--run_qualimap)
    if (run_qualimap) {
        QUALIMAP_BAMQC(bam_bai.map { meta, bam, _bai -> [meta, bam] }, [])
        QUALIMAP_BAMQC_CONSENSUS(consensus_bam_bai.map { meta, bam, _bai -> [meta, bam] }, [])
        ch_multiqc_files = ch_multiqc_files
            .mix(QUALIMAP_BAMQC.out.results.map { _meta, f -> f })
            .mix(QUALIMAP_BAMQC_CONSENSUS.out.results.map { _meta, f -> f })
    }

    ch_multiqc_input = ch_multiqc_files
        .collect()
        .map { files -> [[id: 'multiqc'], files, multiqc_config, [], [], []] }

    MULTIQC(ch_multiqc_input)

    emit:
    report = MULTIQC.out.report // channel: [ val(meta), path(html) ]
    data = MULTIQC.out.data // channel: [ val(meta), path(dir) ]
}
