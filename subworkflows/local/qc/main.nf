//
// QC: FastQC, FastQ Screen, MarkDuplicates metrics, samtools stats/flagstat, Qualimap, the optional
// read info table and the duplex metrics on the mapped and consensus BAMs, summarised in a single MultiQC report.
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
include { GET_READ_INFO } from '../../../modules/local/getreadinfo/main'
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
    output_read_info // value: boolean, write the read info table and compute the duplex metrics from it
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

    // duplex metrics are computed from the mapped BAM (as in Snakemake), not the consensus BAM;
    // with output_read_info the read info table is made first and used as the input instead
    ch_mapped_bam = bam_bai.map { meta, bam, _bai -> [meta, bam] }
    if (output_read_info) {
        GET_READ_INFO(ch_mapped_bam)
        ch_duplex_metrics_input = GET_READ_INFO.out.rinfo
    }
    else {
        ch_duplex_metrics_input = ch_mapped_bam
    }

    CALC_DUPLEX_METRICS(
        ch_duplex_metrics_input,
        fasta_fai.map { meta, fasta, _fai -> [meta, fasta] },
        target_bed,
        output_read_info ? 'rinfo' : 'bam',
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

    // MultiQC exits 0 when it does not know about the duplex plugin (for example a container without
    // duplex-multiqc), and --strict only catches a plugin that crashes. Either way the duplex section
    // would silently be missing, so fail when the plugin wrote none of its files.
    MULTIQC.out.data.subscribe { _meta, data_dir ->
        if (!workflow.stubRun && !data_dir.listFiles().any { f -> f.name.startsWith('duplex_seq') }) {
            error("The MultiQC report has no duplex section: the MULTIQC container must include the duplex-multiqc plugin (see conf/modules.config)")
        }
    }

    emit:
    report = MULTIQC.out.report // channel: [ val(meta), path(html) ]
    data = MULTIQC.out.data // channel: [ val(meta), path(dir) ]
}
