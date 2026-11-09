# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.2.0] - Unreleased

The Nextflow rewrite of the pipeline. The Snakemake pipeline (release 0.1.0) is no longer in the repository; it stays
available at the tag `0.1.0` (`git clone --branch 0.1.0 https://github.com/WEHIGenomicsRnD/duplex-seq-pipeline.git`) and is
the legacy version. The rewrite plan is in [`plans/nextflow-rewrite.md`](plans/nextflow-rewrite.md).

### Added

- A Nextflow (DSL2, nf-core style) implementation of the whole pipeline: `main.nf`, `workflows/duplexseq.nf`,
  `nextflow.config`, `conf/`, nf-core modules (`modules.json`) and local modules.
  - Subworkflows: `PREPARE_REFERENCE` (faidx, sequence dictionary, bwa index; indexes that are supplied are not rebuilt;
    a gzipped FASTA is decompressed first), `FASTQ_ALIGN_UMI` (FastqToSam, fgbio ExtractUmisFromBam, optional
    CorrectUmis), `BAM_ALIGN` (SamToFastq, bwa mem, MergeBamAlignment, index; also used for the consensus reads),
    `CONSENSUS`, `QC` and `VARIANT_CALLING`.
  - Input is a samplesheet (`sample,fastq_1,fastq_2` and the optional `read_structure` and `umi_file` columns), validated
    with nf-schema together with all parameters (`nextflow_schema.json`, `assets/schema_input.json`). Relative paths are
    resolved against the folder of the samplesheet.
  - Profiles: `test`, `docker`, `singularity`, `apptainer`, `conda`, `mamba`, `wave`, `debug`. Tool resources follow
    `config/cluster.yaml` of 0.1.0 (`conf/base.config`), and the hard-coded tool arguments of the Snakemake rules are set in
    `conf/modules.config`, where they can be overridden.
- Variant calling with fgbio ClipBam, nvc and `bin/filter_vcf.py` (local modules `FGBIO_CLIPBAM`, `NVC`, `FILTER_VCF`).
- `--output_read_info` and the `GET_READ_INFO` module; `CALC_DUPLEX_METRICS` reads either the BAM or the read info table.
  The two give identical metrics on the test data.
- A single MultiQC report (`QC/multiqc`) that replaces the two reports of Snakemake. It uses the
  [duplex-multiqc](https://github.com/WEHIGenomicsRnD/duplex-multiqc) plugin 0.1.1 on MultiQC 1.35 (container
  `ghcr.io/wehigenomicsrnd/duplex-multiqc`) and the packaged `assets/multiqc_config.yml`, which lists the duplex metrics
  first and folds each sample's rows into one expandable row in the General Statistics table: the sample row shows the
  duplex metrics, and the fold-out shows the FastQC rows (R1, R2), the raw BAM and the consensus BAM. This needs the stats
  outputs to be named `<sample>_reads_{1,2}` (FastQC), `<sample>_mapped` (samtools stats and flagstat, Qualimap) and
  `<sample>_consensus`. `--multiqc_config` is applied on top of it.
- Guards against a silently missing duplex section in the MultiQC report (a container without the plugin exits 0):
  MultiQC runs with `--strict`, and the pipeline fails if the plugin wrote none of its files.
- `--fastq_screen_db`, a FastQ Screen database folder (see below), and `tests/scripts/build_fastq_screen_db.sh`, which
  builds a one-genome test database from a FASTA file with `bowtie2-build`. The database is not part of the repository.
- Test data (`tests/data`: a 1 kb window of E. coli ATCC 10798 NanoSeq reads, the GCA_002087985.1 reference and a target
  BED), the `test` profile and nf-test tests (`nf-test.config`, `tests/default.nf.test`), including one for FastQ Screen
  that needs `FASTQ_SCREEN_DB`.
- A GitHub Actions workflow (`.github/workflows/nf-test.yml`) that runs the nf-tests with docker on pushes and pull requests
  to `main` and `dev` only. It first builds the FastQ Screen test database with bowtie2 from `tests/scripts/environment.yml`.
- `--merge_lanes`: the rows of a sample in the samplesheet are its lanes, and their files are concatenated in lane order
  into one pair of FASTQ files per sample (nf-core `CAT_FASTQ`, in the subworkflow `MERGE_LANES`). The lane is read from file
  names such as `_L001_R1.fastq.gz`, and the rows are checked before anything runs. Without the option a repeated sample
  name is still an error, but the check is now made at start-up and no longer by the samplesheet schema.
- `--markdup_duplex_umi` (on by default): MarkDuplicates marks duplicates per duplex molecule (`--DUPLEX_UMI true`); turn it
  off for single-UMI read structures or to reproduce the per-strand marking of 0.1.0.
- `meta.yml` for the local modules, in the nf-core format (checked against the nf-core module meta schema).
- Tooling: pre-commit hooks (prettier, editorconfig-checker, `nextflow lint`, schema check) and a lint workflow.

### Changed

Compared with the Snakemake pipeline (0.1.0):

- **Tool versions**, as installed by the Nextflow modules (Snakemake in brackets): fgbio 3.1.2 (1.4.0), Picard 3.5.0 and
  GATK 4.7.0.0 (Picard 2.27.4), bwa 0.7.19 (0.7.17), samtools 1.22.1 and 1.25 (1.14), FastQC 0.13.0 (0.11.8),
  FastQ Screen 0.16.0 (0.14.0), Qualimap 2.3 (2.2.2), MultiQC 1.35 (1.14), nvc 0.0.4 (same). Output can differ slightly
  because of these upgrades.
- **Read group:** FastqToSam uses the sample name as the read group name. Snakemake used `null`, which Picard 2.27 let
  through but which writes an `@RG` header line without tagging any read; Picard 3.5 MarkDuplicates fails on it.
- **Realignment** uses GATK4 `SamToFastq` and `MergeBamAlignment` (there is no nf-core Picard module for the latter), so
  bwa mem reads two FASTQ files instead of one interleaved file, and the BAM index comes from `samtools index`.
- **Consensus calling** follows the fgbio 3 command line: the BAMs are queryname sorted for `MergeBamAlignment`,
  `FilterConsensusReads` and `ClipBam` (fgbio >= 2 requires query-grouped input), `--max-no-call-fraction` is now
  `--max-no-calls`, and the filtered BAM is written in coordinate order.
- **QC:** duplicates are marked per duplex molecule (`--markdup_duplex_umi`, on by default; Snakemake marked per strand,
  see the notes below), samtools stats, flagstat and Qualimap read the duplicate-marked BAM, only the MarkDuplicates
  metrics are published (not the marked BAM), and the duplex metrics are named after the sample (`--sample`), so they
  merge into the sample's row in MultiQC.
- **fgbio memory:** the nf-core fgbio modules cap the Java heap at 8 GB however much memory the task has (they only lower
  it). `GroupReadsByUmi`, `CallDuplexConsensusReads` and `FilterConsensusReads` are patched (`modules.json`) to use 90% of
  the task memory, as the local fgbio modules do, so the 64 GB reserved in `conf/base.config` is used as in 0.1.0.
- **Parameters** keep their Snakemake names and defaults, except `ref` (now `--fasta`), `FastqScreen_config` (now
  `--fastq_screen_db`, a database folder in the `FastQ_Screen_Genomes` layout of `fastq_screen --get_genomes`, because
  the nf-core module stages a whole folder) and `cluster_config` (replaced by Nextflow profiles and resources). New
  parameters are `--input` (the samplesheet), `--outdir`, `--bwa`, `--dict` and `--fasta_fai` (optional pre-built
  indexes), `--multiqc_title`, and the memory and buffer settings that were in `config/cluster.yaml` of 0.1.0
  (`--sort_mem`, `--bam_chunk_size`, `--read_info_buffer_size`, `--markdup_max_records_in_ram`).
- `bin/filter_vcf.py` was rewritten (a command line, `filter_vcf.py INPUT_VCF OUTPUT_VCF`, instead of the Snakemake
  globals; readable names and comments). It gives the same output as the Snakemake script on nvc VCFs, with two
  deliberate differences: other INFO fields are kept when an N allele is removed (the script dropped everything except
  AC and AF), and the DP header line is also added when a VCF has no `##FORMAT` lines.
- The editorconfig check covers the whole repository (except the nf-core directories) and does not check the indent size
  of markdown, because prettier aligns content under numbered list items with 3 spaces. Trailing whitespace in
  the Snakemake `config/config.yaml` and the docstring indentation in `slurm/slurm_utils.py` were fixed (before those files were removed).

### Removed

- varscan is not implemented in the Nextflow pipeline (only nvc), and with it the parameters `variant_caller`,
  `min_coverage` and `min_vaf`. The Snakemake pipeline (0.1.0) still supports it.
- The Snakemake pipeline itself (`workflow/`, `config/`, `slurm/` and `dag.pdf`), and the Snakemake half of the README. The
  tag `0.1.0` has all of it, including its README.
- The unpublished first skeleton (`nf-pipeline/`).

### Notes

#### Duplicate marking and the QC statistics

In 0.1.0 MarkDuplicates ran next to the statistics tools, which read the BAM before the duplicates were marked. In the
rewrite they read the marked BAM. On the 1 kb test window (87.6% duplicates, so the sizes are not typical; the BAMs differ
only by the duplicate flag):

| Tool     | Unmarked BAM (0.1.0)                                  | Marked per strand (`--markdup_duplex_umi false`) | Marked per molecule (default)          |
| -------- | ----------------------------------------------------- | ------------------------------------------------ | -------------------------------------- |
| flagstat | 0 duplicates                                          | 71,041 (82.29% of primary mapped reads)          | 75,637 (87.61%)                        |
| stats    | `reads duplicated` 0; MAPQ histogram of all 86,329    | 71,041; MAPQ histogram of 15,288 reads           | 75,637; MAPQ histogram of 10,692 reads |
| Qualimap | its own position-only estimate, 84,796 reads (96.85%) | the flagged duplicates, 71,041                   | the flagged duplicates, 75,637         |

- The flagstat count is exactly the Picard figure (per strand: 2 x 35,443 pair duplicates + 155 unpaired duplicates,
  `PERCENT_DUPLICATION` 0.82291; per molecule: 2 x 37,747 + 143, 0.876148). Coverage, insert size, GC, read length, error
  rate and mapping quality are the same whichever BAM is used. The MAPQ histogram of samtools stats leaves out reads
  that are flagged as duplicates.
- Qualimap cannot see UMIs, so on the unmarked BAM its duplication rate contradicts Picard's.
- Nothing downstream is affected by the flags: `calc-duplex-metrics` gives identical output on both BAMs,
  `fgbio GroupReadsByUmi` and `get_read_info.py` do not filter on the flag. Even so, consensus calling and the duplex
  metrics keep the unmarked BAM, so that no tool that honours the flag can drop members of a UMI family.
- The statistics now wait for MarkDuplicates, and a MarkDuplicates failure also stops them. The marked BAM was already
  written (and discarded) before, so there is no extra computation.
- A parity run against Snakemake shows these duplicate lines (flagstat, stats, Qualimap, Picard) differing; they match
  when the rewrite is run with `--markdup_duplex_umi false`, apart from the flagstat, stats and Qualimap lines that now
  count the flags.

#### Duplex UMIs in MarkDuplicates

MarkDuplicates is run with `--BARCODE_TAG RX`. The `RX` tag is `A-B` on one strand of a molecule and `B-A` on the other,
so by default (and in 0.1.0) the two strands count as different molecules and the flags mark duplicates per strand.
Picard 3.5 can treat them as one molecule with `--DUPLEX_UMI true`, which the rewrite now does by default. On the same
test BAM:

| MarkDuplicates                | Pair duplicates | % duplication | Unique read pairs | Library size estimate |
| ----------------------------- | --------------- | ------------- | ----------------- | --------------------- |
| per strand (0.1.0)            | 35,443          | 82.3%         | 7,600             | 7,628                 |
| per molecule (`--DUPLEX_UMI`) | 37,747          | 87.6%         | 5,296             | 5,297                 |

The unique read pairs are close to the family counts of `calc-duplex-metrics` on the same BAM: 7,880 strand-level
families (2 x 2,315 paired + 3,250 single-strand) against 7,600 per strand, and 5,565 molecule-level families
(`total_families`) against 5,296 per molecule. So the per-molecule duplicates are counted in the unit of the duplex
metrics.

`--DUPLEX_UMI` needs every `RX` tag to be two equal-length UMIs joined by a hyphen: Picard stops with "Duplex UMIs must be
of the form X-Y where X and Y are equal length UMIs" otherwise (a missing `RX` tag is not an error). Use
`--markdup_duplex_umi false` for read structures with a single UMI, and to reproduce the per-strand marking of 0.1.0.

#### Parity with 0.1.0

Both pipelines were run on the 1 kb test window (sample `NanoMB1Rep1`): Snakemake 8.11.3 with its own conda environments
and this branch with apptainer, with `--markdup_duplex_umi false` and with the default (set-up and comparison script in
`tests/parity`). The Nextflow pipeline gives the same results in the early steps and differs from fgbio 3 onwards.

**Identical:**

- the aligned BAM, in all 86,593 records (position, MAPQ, CIGAR, mate fields, sequence, qualities and the UMI tags), so
  bwa 0.7.17 and 0.7.19 align identically. Only the `RG` tag differs (deliberate, see above) and the `MQ` tag, which is
  missing on the 21 supplementary records in 0.1.0 (not isolated; probably GATK 4.7 against Picard 2.27 MergeBamAlignment);
- the FastQ Screen results, the FastQC basic statistics, the (empty) filtered VCF and 19 of the 20 duplex metrics (to
  1e-9);
- the Picard MarkDuplicates metrics per strand (35,443 pair duplicates, 155 unpaired, 0.82291) apart from the optical
  duplicates.

**Different, with the cause** (each tested by running the tools of both versions on the same input):

- **UMI grouping (fgbio 1.4.0 against 3.1.2):** 85,572 records and 4,254 molecule families against 86,304 records and
  4,359 families. fgbio 3 also groups the 366 templates (732 reads) that have an unmapped read (241) or mates on different
  contigs (125); fgbio 1.4.0 drops them. All 4,254 families of 0.1.0 are unchanged.
- **Consensus calling (fgbio 1.4.0 against 3.1.2):** 4,690 against 4,694 consensus reads. The new grouping adds 8; the
  consensus caller itself gives a different sequence for about 16% of the reads (3,947 of 4,690 sequences are identical)
  and different qualities for more, with the same mean base quality (87.8 against 87.9). Each pipeline's consensus was
  reproduced exactly with the standalone fgbio (4,690 and 4,694 reads).
- **Filtered consensus reads:** 3,380 against 3,580. `FilterConsensusReads` gives the same result in both versions on the
  same input, so this follows from the consensus calls.
- **nvc:** 801 against 761 sites (697 shared). Every site in both is an N-only call, so no variant is called and both
  filtered VCFs are empty. The sites that differ are consensus bases that one version calls as N.
- **Duplex metrics:** Snakemake has one more metric, `unpaired_families`. It installs `calculate-duplex-metrics` from the
  repository at run time (commit `39fc6c9`), which has the metric; the `0.2.0` container used in this run did not. The pipeline
  now uses `calculate-duplex-metrics` 0.2.1, which has it: on the same BAM all 20 metrics are identical to the Snakemake run
  (checked directly, and the nf-tests give 20 metrics; the comparison itself was not repeated).
- **samtools stats:** 101 pairs that samtools 1.14 reports as outward oriented are inward oriented in 1.23 and later (the
  same on either BAM); the duplicate counts differ because of the marked BAM (see above).
- **Picard:** 388 against 390 optical duplicates, from Picard 2.27.4 against 3.5.0 (2.27.4 gives 388 on the Nextflow BAM
  too). With the default `--markdup_duplex_umi` the duplicate figures differ further, as described above.
- **Qualimap:** the duplication lines (flagged against estimated), two extra lines for overlapping read pairs (the nf-core
  module runs Qualimap with `--collect-overlap-pairs`), and for the consensus BAM the consequences of the consensus
  differences.

Decision: fgbio 3.1.2 is kept (as in the nf-core modules), to take improvements where they come. The consensus, and with
it the filtered BAM and the nvc calls, therefore do not match 0.1.0 exactly; exact parity would need fgbio 1.4.0 in the
local and nf-core fgbio modules. The duplex metrics container is `calculate-duplex-metrics` 0.2.1, which has `unpaired_families`.

## [0.1.0]

The Snakemake pipeline (tag `0.1.0`, commit `d71e733`): FASTQ files to UMI-extracted and aligned BAM, duplex consensus
calling, QC (FastQC, FastQ Screen, MarkDuplicates, samtools, Qualimap, `calc-duplex-metrics` and MultiQC with the
duplex-multiqc plugin) and variant calling (nvc or varscan).
