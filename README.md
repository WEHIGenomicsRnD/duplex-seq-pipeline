# Duplex Sequencing Pipeline

A [Nextflow](https://www.nextflow.io) pipeline for processing duplex
sequencing data (e.g.
[NanoSeq](https://www.nature.com/articles/s41586-021-03477-4)) to produce
duplex consensus reads, variant calls, and QC metrics. The workflow follows the
IDT recommended approach for [processing sequence data with unique molecular
identifiers](https://plone.bcgsc.ca/services/solseq/duplex-umi-documents/idt_analysisguideline_varcall-umis-dupseqadapters/).

This is the Nextflow rewrite of the pipeline (release 0.2.0, not released yet). It runs the same analysis as the original
[Snakemake version](#legacy-snakemake-version) (release 0.1.0) with containers or conda. On the test data both give the
same alignments, and the remaining differences were traced to tool versions, mainly fgbio 3 against fgbio 1.4.0: see
[`tests/parity`](tests/parity/README.md) and the [changelog](CHANGELOG.md).

## What it does

1. **Reference:** a samtools index, a sequence dictionary and a bwa index are built unless you supply them. A gzipped FASTA is
   decompressed first.
2. **UMIs and alignment:** the FASTQ files become an unmapped BAM, the UMIs are extracted with fgbio (`--read_structure`) and
   optionally corrected against a list of known UMIs (`--umis`), and the reads are aligned with bwa mem and merged back.
3. **Duplex consensus:** reads are grouped by UMI, duplex consensus reads are called, realigned and filtered (fgbio).
4. **Variants:** overlapping read pairs of the filtered consensus reads are clipped (fgbio ClipBam), variants are called with
   the naive variant caller (nvc) and the VCF is filtered.
5. **QC:** FastQC, FastQ Screen (when a database is given), Picard MarkDuplicates (per duplex molecule), samtools stats and
   flagstat, Qualimap (`--run_qualimap`), the coverage of the duplex reads (mosdepth, on the clipped consensus BAM) and the
   duplex metrics of
   [calc-duplex-metrics](https://github.com/WEHIGenomicsRnD/calculate-duplex-metrics), all in a single MultiQC report that
   includes the duplex metrics (via the [duplex-multiqc](https://github.com/WEHIGenomicsRnD/duplex-multiqc) plugin).

## Running

```bash
# check the installation on the small test data set
nextflow run . -profile test,docker --outdir results

# run on your data
nextflow run . -profile docker --input samplesheet.csv --fasta genome.fasta --outdir results
```

You need Nextflow 25.10 or later and a container engine: use the profile `docker`, `singularity` or `apptainer` (`wave`
builds containers on demand). The `conda` and `mamba` profiles cannot run the duplex metrics, because
calc-duplex-metrics is not on conda and only comes as a container
(`ghcr.io/wehigenomicsrnd/calculate-duplex-metrics`). Run `nextflow run . --help` for all parameters.

The samplesheet is a CSV with the columns below. Leave a cell empty to use the default; do not write `null` or `NA`, which
would be read as a file name. Relative paths are resolved against the folder of the samplesheet.

| Column           | Required | Description                                                                                                |
| ---------------- | -------- | ---------------------------------------------------------------------------------------------------------- |
| `sample`         | yes      | Sample name, without spaces.                                                                               |
| `fastq_1`        | yes      | FASTQ file of read 1 (`.fastq.gz` or `.fq.gz`).                                                            |
| `fastq_2`        | yes      | FASTQ file of read 2.                                                                                      |
| `read_structure` | no       | fgbio read structure of this sample, one per read (e.g. `3M2S146T 3M2S146T`). Default: `--read_structure`. |
| `umi_file`       | no       | File of known UMIs; the UMIs of this sample are corrected against it. Default: `--umis`, or no correction. |

```csv
sample,fastq_1,fastq_2
NanoMB1Rep1,NanoMB1Rep1_R1.fastq.gz,NanoMB1Rep1_R2.fastq.gz
```

Read structures and UMI correction are decided per sample, so one run can mix kits. In the example below the first sample is
corrected against `umis.txt` and the second is not corrected and has another read structure:

```csv
sample,fastq_1,fastq_2,read_structure,umi_file
A,A_R1.fastq.gz,A_R2.fastq.gz,3M2S146T 3M2S146T,umis.txt
B,B_R1.fastq.gz,B_R2.fastq.gz,3M2S+T 3M2S+T,
```

A sample is corrected when it has a `umi_file` or, if that cell is empty, when `--umis` is given. An empty cell cannot switch
off a global `--umis`, so for a run in which only some samples are corrected, do not set `--umis` and give the `umi_file` to
those samples.

If a sample was sequenced on several lanes, give every lane its own row with the same sample name and use `--merge_lanes`.
The files of the lanes are then concatenated, in lane order, into one pair of FASTQ files per sample before anything else
runs. The lane is read from the file names, which have to end like `_L001_R1.fastq.gz` and `_L001_R2.fastq.gz` (any lane
number; the bcl2fastq form `_L001_R1_001.fastq.gz` works too):

```csv
sample,fastq_1,fastq_2
NanoMB1Rep1,NanoMB1Rep1_L001_R1.fastq.gz,NanoMB1Rep1_L001_R2.fastq.gz
NanoMB1Rep1,NanoMB1Rep1_L002_R1.fastq.gz,NanoMB1Rep1_L002_R2.fastq.gz
```

Without `--merge_lanes` a sample name on more than one row is an error. With it, the run stops at the start if a lane
cannot be read from a file name, if R1 and R2 of a row are from different lanes, if a lane is listed twice for a sample, or
if the rows of a sample have different read structures or UMI files. A sample on a single row is left as it is. The merged
files are only kept in the work directory.

The main parameters (all in `nextflow_schema.json`):

| Parameter                             | Default             | Description                                                                                                                                                                                                                             |
| ------------------------------------- | ------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--merge_lanes`                       | `false`             | The rows of a sample in the samplesheet are its lanes: merge them (see above).                                                                                                                                                          |
| `--input`, `--outdir`                 |                     | The samplesheet and the output directory (`results`).                                                                                                                                                                                   |
| `--fasta`                             |                     | Reference genome FASTA (may be gzipped). `--fasta_fai`, `--dict` and `--bwa` take indexes that already exist.                                                                                                                           |
| `--read_structure`                    | `3M2S146T 3M2S146T` | fgbio read structure of the two reads.                                                                                                                                                                                                  |
| `--umis`                              |                     | Text file of known UMIs (one per line); UMIs are corrected against it. Without it no correction is done.                                                                                                                                |
| `--target_bed`                        |                     | BED file of target regions; adds the on-target rates to the duplex metrics.                                                                                                                                                             |
| `--fastq_screen_db`                   |                     | FastQ Screen database folder: the `FastQ_Screen_Genomes` layout of `fastq_screen --get_genomes`, with a `fastq_screen.conf`. FastQ Screen is skipped without it. `tests/scripts/build_fastq_screen_db.sh` builds a small one for tests. |
| `--run_qualimap`                      | `false`             | Run Qualimap (resource-intensive).                                                                                                                                                                                                      |
| `--output_read_info`                  | `false`             | Write the read info table (`QC/read_info`) and compute the duplex metrics from it instead of from the BAM.                                                                                                                              |
| `--markdup_duplex_umi`                | `true`              | Mark duplicates per duplex molecule (the UMIs `A-B` and `B-A` of the two strands count as one). Turn it off for single-UMI read structures, or to mark per strand as release 0.1.0 does.                                                |
| `--multiqc_config`, `--multiqc_title` |                     | A MultiQC config that is applied on top of the packaged one, and a report title.                                                                                                                                                        |

The consensus calling and filtering options (`--min_reads`, `--min_base_quality`, ...) and the variant calling options
(`--min_reads2`, `--ploidy`) keep the names and defaults of `config/config.yaml`; `ref` is now `--fasta`, and
`FastqScreen_config` is now `--fastq_screen_db`. varscan is not implemented, only nvc.

### FastQ Screen database

`--fastq_screen_db` is a folder holding a `fastq_screen.conf` and the bowtie2 indexes it lists. The easiest way to get a
full set is `fastq_screen --get_genomes`, which downloads the genomes and writes the config. To screen against your own
genomes, write the config yourself: the
[example config](https://github.com/StevenWingett/FastQ-Screen/blob/master/fastq_screen.conf.example) from the FastQ Screen
repository shows the `DATABASE` lines. Give the index paths relative to the folder (for example
`DATABASE	Human	FastQ_Screen_Genomes/Human/Homo_sapiens.GRCh38`) or as absolute paths that are visible inside the container. The
config of the folder is the one used, so keep one folder per set of genomes.

## Outputs

| Directory         | Content                                                                                                                                                                                                                           |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `reference/`      | The indexes that were built (`.fai`, `.dict`, `bwa/`).                                                                                                                                                                            |
| `mapped_bams/`    | `<sample>_mapped_merged.bam`: the aligned BAM with the UMI tags.                                                                                                                                                                  |
| `grouped_by_umi/` | The BAM grouped by UMI and the family size histogram.                                                                                                                                                                             |
| `consensus/`      | The duplex consensus reads: unmapped, mapped, filtered (`_filtered.bam`) and clipped for variant calling (`_filtered_clipped.bam`).                                                                                               |
| `variants/`       | `<sample>.vcf` from nvc and `<sample>_filtered.vcf`.                                                                                                                                                                              |
| `QC/`             | FastQC, FastQ Screen, MarkDuplicates metrics, samtools, Qualimap, `mosdepth/` (coverage of the duplex reads), `duplex_metrics/` and, with `--output_read_info`, `read_info/`. The consensus reads have their own `QC/consensus/`. |
| `QC/multiqc/`     | The MultiQC report (`multiqc_report.html`).                                                                                                                                                                                       |
| `pipeline_info/`  | Tool versions (`software_versions.yml`), the execution report, timeline, trace and DAG.                                                                                                                                           |

## Development

The tests use [nf-test](https://www.nf-test.com): `nf-test test --profile test,docker`. The FastQ Screen test needs a
database, so run it with `FASTQ_SCREEN_DB` set to the `FastQ_Screen_Genomes` folder made by
`tests/scripts/build_fastq_screen_db.sh` (bowtie2 is in `tests/scripts/environment.yml`). A GitHub Actions workflow
(`.github/workflows/nf-test.yml`) builds that database and runs all tests with docker on pushes and pull requests to `main`
and `dev`. The linters (prettier, editorconfig,
`nextflow lint`) run as [pre-commit](https://pre-commit.com) hooks. `plans/nextflow-rewrite.md` is the plan of the rewrite and
`CHANGELOG.md` lists the changes with respect to the Snakemake version (0.1.0).

## Legacy Snakemake version

The first release of the pipeline (0.1.0) was written in [Snakemake](https://snakemake.readthedocs.io). It is no longer
maintained or part of this branch, but you can still run it from the tag, which has its own README with the
installation, configuration (`config/config.yaml`, `config/fastq_screen.conf` of that tag) and output description:

```bash
git clone --branch 0.1.0 https://github.com/WEHIGenomicsRnD/duplex-seq-pipeline.git duplex-seq-pipeline-0.1.0
cd duplex-seq-pipeline-0.1.0
mamba create -c conda-forge -c bioconda -n snakemake snakemake=7.0.0-0
conda activate snakemake
snakemake --use-conda --conda-frontend mamba --cores 4
```

The Nextflow pipeline gives the same alignments; the differences in the other results and the changed defaults are listed
in the [changelog](CHANGELOG.md), and [`tests/parity`](tests/parity/README.md) compares the two on the test data.
