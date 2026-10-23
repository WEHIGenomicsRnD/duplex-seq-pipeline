#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    duplex-seq-pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Duplex sequencing pipeline for consensus calling, QC and variant detection.
    Converted from Snakemake to Nextflow DSL2 with nf-core modular structure.

    https://github.com/your-org/duplex-seq-pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// -- nf-core modules (reusable from the community registry) --
include { PICARD_CREATESEQUENCEDICTIONARY            } from './modules/nf-core/picard/createsequencedictionary/main'
include { BWA_INDEX                                   } from './modules/nf-core/bwa/index/main'
include { BWA_MEM                                     } from './modules/nf-core/bwa/mem/main'
include { BWA_MEM as BWA_MEM_CONSENSUS                } from './modules/nf-core/bwa/mem/main'
include { FASTQC                                      } from './modules/nf-core/fastqc/main'
include { SAMTOOLS_STATS                              } from './modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT                           } from './modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_CONSENSUS  } from './modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_CONSENSUS } from './modules/nf-core/samtools/flagstat/main'
include { PICARD_MARKDUPLICATES                       } from './modules/nf-core/picard/markduplicates/main'
include { QUALIMAP_BAMQC                              } from './modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC as QUALIMAP_BAMQC_CONSENSUS  } from './modules/nf-core/qualimap/bamqc/main'
include { MULTIQC                                     } from './modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_CONSENSUS                } from './modules/nf-core/multiqc/main'
include { FGBIO_GROUPREADSBYUMI                       } from './modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_CALLDUPLEXCONSENSUSREADS              } from './modules/nf-core/fgbio/callduplexconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS                  } from './modules/nf-core/fgbio/filterconsensusreads/main'

// -- local modules --
include { PICARD_FASTQTOSAM                           } from './modules/local/picard/fastqtosam/main'
include { FGBIO_EXTRACT_UMIS                          } from './modules/local/fgbio/extract_umis/main'
include { FGBIO_CORRECT_UMIS                          } from './modules/local/fgbio/correct_umis/main'
include { PICARD_SAMTOFASTQ                           } from './modules/local/picard/samtofastq/main'
include { PICARD_SAMTOFASTQ as PICARD_SAMTOFASTQ_CONSENSUS } from './modules/local/picard/samtofastq/main'
include { PICARD_MERGEBAMALIGNMENT                    } from './modules/local/picard/mergebamalignment/main'
include { PICARD_MERGEBAMALIGNMENT as PICARD_MERGEBAMALIGNMENT_CONSENSUS } from './modules/local/picard/mergebamalignment/main'
include { FASTQ_SCREEN                                } from './modules/local/fastq_screen/main'
include { GET_READ_INFO                               } from './modules/local/get_read_info/main'
include { FORMAT_READ_INFO                            } from './modules/local/format_read_info/main'
include { CALC_DUPLEX_METRICS                         } from './modules/local/calc_duplex_metrics/main'

// -- local subworkflows --
include { VARIANT_CALLING                             } from './subworkflows/local/variant_calling'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // -- Validate inputs --
    if (!params.input) {
        error "Please provide a samplesheet via --input"
    }
    if (!params.ref) {
        error "Please provide a reference genome via --ref"
    }

    def ch_versions = Channel.empty()

    // -----------------------------------------------------------------------
    //  Parse samplesheet
    // -----------------------------------------------------------------------
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            [meta, [file(row.fastq_1), file(row.fastq_2)]]
        }
        .set { ch_reads }

    // -----------------------------------------------------------------------
    //  Prepare reference
    // -----------------------------------------------------------------------
    def ch_ref  = Channel.fromPath(params.ref, checkIfExists: true).collect()
    def ch_fai  = ch_ref   // assumed to exist alongside .fasta

    PICARD_CREATESEQUENCEDICTIONARY(ch_ref.map { ref -> [[id: 'genome'], ref] })
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)
    def ch_dict = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict.map { meta, dict -> dict }.collect()

    BWA_INDEX(ch_ref.map { ref -> [[id: 'genome'], ref] })
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    def ch_bwa_index = BWA_INDEX.out.index

    // -----------------------------------------------------------------------
    //  PREPROCESSING: FASTQ → unmapped BAM → UMI extraction → align → merge
    // -----------------------------------------------------------------------
    def correct_umis = params.umis != ''

    // Convert FASTQ to unmapped BAM
    PICARD_FASTQTOSAM(ch_reads)
    ch_versions = ch_versions.mix(PICARD_FASTQTOSAM.out.versions.first())

    // Extract UMIs
    FGBIO_EXTRACT_UMIS(PICARD_FASTQTOSAM.out.bam)
    ch_versions = ch_versions.mix(FGBIO_EXTRACT_UMIS.out.versions.first())

    // Optionally correct UMIs
    if (correct_umis) {
        def ch_umis = Channel.fromPath(params.umis, checkIfExists: true).collect()
        FGBIO_CORRECT_UMIS(FGBIO_EXTRACT_UMIS.out.bam, ch_umis)
        ch_versions = ch_versions.mix(FGBIO_CORRECT_UMIS.out.versions.first())
        ch_umi_bam = FGBIO_CORRECT_UMIS.out.bam
    } else {
        ch_umi_bam = FGBIO_EXTRACT_UMIS.out.bam
    }

    // BAM → interleaved FASTQ
    PICARD_SAMTOFASTQ(ch_umi_bam)
    ch_versions = ch_versions.mix(PICARD_SAMTOFASTQ.out.versions.first())

    // Align with BWA MEM (interleaved input via ext.args = '-p')
    BWA_MEM(PICARD_SAMTOFASTQ.out.fastq, ch_bwa_index, ch_ref, true)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    // Merge aligned + unmapped BAMs to restore UMI tags
    def ch_merge_input = ch_umi_bam
        .join(BWA_MEM.out.bam)
        .map { meta, unmapped, mapped -> [meta, unmapped, mapped] }

    PICARD_MERGEBAMALIGNMENT(ch_merge_input, ch_ref, ch_dict)
    ch_versions = ch_versions.mix(PICARD_MERGEBAMALIGNMENT.out.versions.first())

    // -----------------------------------------------------------------------
    //  QC — raw reads
    // -----------------------------------------------------------------------
    FASTQC(ch_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    if (params.fastq_screen_config) {
        def ch_fqs_config = Channel.fromPath(params.fastq_screen_config, checkIfExists: true).collect()
        FASTQ_SCREEN(ch_reads, ch_fqs_config)
        ch_versions = ch_versions.mix(FASTQ_SCREEN.out.versions.first())
    }

    PICARD_MARKDUPLICATES(PICARD_MERGEBAMALIGNMENT.out.bam, ch_ref, ch_dict)
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    SAMTOOLS_STATS(PICARD_MERGEBAMALIGNMENT.out.bam, ch_ref)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    SAMTOOLS_FLAGSTAT(PICARD_MERGEBAMALIGNMENT.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    QUALIMAP_BAMQC(PICARD_MERGEBAMALIGNMENT.out.bam, [])
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

    // Read info for duplex metrics (on raw merged BAM)
    GET_READ_INFO(PICARD_MERGEBAMALIGNMENT.out.bam)
    ch_versions = ch_versions.mix(GET_READ_INFO.out.versions.first())

    FORMAT_READ_INFO(GET_READ_INFO.out.rinfo)
    ch_versions = ch_versions.mix(FORMAT_READ_INFO.out.versions.first())

    CALC_DUPLEX_METRICS(FORMAT_READ_INFO.out.formatted, ch_ref)
    ch_versions = ch_versions.mix(CALC_DUPLEX_METRICS.out.versions.first())

    // MultiQC — raw reads
    def ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { meta, zip -> zip })
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect { meta, m -> m })
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect { meta, s -> s })
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect { meta, f -> f })
    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect { meta, r -> r })

    def ch_multiqc_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true).collect() :
        Channel.empty()

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config,
        Channel.empty(),
        Channel.empty()
    )
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    // -----------------------------------------------------------------------
    //  CONSENSUS: group by UMI → duplex consensus → remap → filter
    // -----------------------------------------------------------------------
    FGBIO_GROUPREADSBYUMI(PICARD_MERGEBAMALIGNMENT.out.bam)
    ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions.first())

    FGBIO_CALLDUPLEXCONSENSUSREADS(FGBIO_GROUPREADSBYUMI.out.bam)
    ch_versions = ch_versions.mix(FGBIO_CALLDUPLEXCONSENSUSREADS.out.versions.first())

    // Remap consensus reads
    PICARD_SAMTOFASTQ_CONSENSUS(FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam)
    ch_versions = ch_versions.mix(PICARD_SAMTOFASTQ_CONSENSUS.out.versions.first())

    BWA_MEM_CONSENSUS(PICARD_SAMTOFASTQ_CONSENSUS.out.fastq, ch_bwa_index, ch_ref, true)
    ch_versions = ch_versions.mix(BWA_MEM_CONSENSUS.out.versions.first())

    // Merge consensus aligned + unmapped
    def ch_consensus_merge = FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam
        .join(BWA_MEM_CONSENSUS.out.bam)
        .map { meta, unmapped, mapped -> [meta, unmapped, mapped] }

    PICARD_MERGEBAMALIGNMENT_CONSENSUS(ch_consensus_merge, ch_ref, ch_dict)
    ch_versions = ch_versions.mix(PICARD_MERGEBAMALIGNMENT_CONSENSUS.out.versions.first())

    // Filter consensus reads
    FGBIO_FILTERCONSENSUSREADS(PICARD_MERGEBAMALIGNMENT_CONSENSUS.out.bam, ch_ref)
    ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions.first())

    // -----------------------------------------------------------------------
    //  QC — consensus reads
    // -----------------------------------------------------------------------
    SAMTOOLS_STATS_CONSENSUS(FGBIO_FILTERCONSENSUSREADS.out.bam, ch_ref)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_CONSENSUS.out.versions.first())

    SAMTOOLS_FLAGSTAT_CONSENSUS(FGBIO_FILTERCONSENSUSREADS.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_CONSENSUS.out.versions.first())

    QUALIMAP_BAMQC_CONSENSUS(FGBIO_FILTERCONSENSUSREADS.out.bam, [])
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC_CONSENSUS.out.versions.first())

    // MultiQC — consensus reads
    def ch_multiqc_consensus_files = Channel.empty()
    ch_multiqc_consensus_files = ch_multiqc_consensus_files.mix(
        SAMTOOLS_STATS_CONSENSUS.out.stats.collect { meta, s -> s }
    )
    ch_multiqc_consensus_files = ch_multiqc_consensus_files.mix(
        SAMTOOLS_FLAGSTAT_CONSENSUS.out.flagstat.collect { meta, f -> f }
    )
    ch_multiqc_consensus_files = ch_multiqc_consensus_files.mix(
        QUALIMAP_BAMQC_CONSENSUS.out.results.collect { meta, r -> r }
    )

    MULTIQC_CONSENSUS(
        ch_multiqc_consensus_files.collect(),
        ch_multiqc_config,
        Channel.empty(),
        Channel.empty()
    )

    // -----------------------------------------------------------------------
    //  VARIANT CALLING
    // -----------------------------------------------------------------------
    VARIANT_CALLING(
        FGBIO_FILTERCONSENSUSREADS.out.bam,
        ch_ref,
        ch_dict,
        params.variant_caller
    )
    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
