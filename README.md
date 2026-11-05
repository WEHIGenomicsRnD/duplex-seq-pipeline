# Duplex Sequencing Pipeline

A [Snakemake](https://snakemake.readthedocs.io) pipeline for processing duplex
sequencing data (e.g.
[NanoSeq](https://www.nature.com/articles/s41586-021-03477-4)) to produce
duplex consensus reads, variant calls, and QC metrics. The workflow follows the
IDT recommended approach for [processing sequence data with unique molecular
identifiers](https://plone.bcgsc.ca/services/solseq/duplex-umi-documents/idt_analysisguideline_varcall-umis-dupseqadapters/).

## Nextflow version (in development)

A Nextflow rewrite of the pipeline (release 0.2.0, not released yet) runs the same analysis with containers or conda. The
[Snakemake version](#snakemake-version) (release 0.1.0, described further down) remains the reference implementation.
On the test data both give the same alignments, and the remaining differences were traced to tool versions, mainly fgbio 3
against fgbio 1.4.0: see [`tests/parity`](tests/parity/README.md) and the [changelog](CHANGELOG.md).

### What it does

1. **Reference:** a samtools index, a sequence dictionary and a bwa index are built unless you supply them. A gzipped FASTA is
   decompressed first.
2. **UMIs and alignment:** the FASTQ files become an unmapped BAM, the UMIs are extracted with fgbio (`--read_structure`) and
   optionally corrected against a list of known UMIs (`--umis`), and the reads are aligned with bwa mem and merged back.
3. **Duplex consensus:** reads are grouped by UMI, duplex consensus reads are called, realigned and filtered (fgbio).
4. **Variants:** overlapping read pairs of the filtered consensus reads are clipped (fgbio ClipBam), variants are called with
   the naive variant caller (nvc) and the VCF is filtered.
5. **QC:** FastQC, FastQ Screen (when a database is given), Picard MarkDuplicates (per duplex molecule), samtools stats and
   flagstat, Qualimap (`--run_qualimap`) and the duplex metrics of
   [calc-duplex-metrics](https://github.com/WEHIGenomicsRnD/calculate-duplex-metrics), all in a single MultiQC report that
   includes the duplex metrics (via the [duplex-multiqc](https://github.com/WEHIGenomicsRnD/duplex-multiqc) plugin).

### Running

```bash
# check the installation on the small test data set
nextflow run . -profile test,docker --outdir results

# run on your data
nextflow run . -profile docker --input samplesheet.csv --fasta genome.fasta --outdir results
```

You need Nextflow 25.04 or later and a container engine: use the profile `docker`, `singularity` or `apptainer` (`wave`
builds containers on demand). The `conda` and `mamba` profiles cannot run the duplex metrics, because
calc-duplex-metrics is not on conda and only comes as a container
(`ghcr.io/wehigenomicsrnd/calculate-duplex-metrics`). Run `nextflow run . --help` for all parameters.

The samplesheet is a CSV with columns `sample,fastq_1,fastq_2` and the optional
columns `read_structure` and `umi_file`, which override `--read_structure` and
`--umis` for that sample. Relative paths are resolved against the folder of the samplesheet:

```csv
sample,fastq_1,fastq_2
NanoMB1Rep1,NanoMB1Rep1_R1.fastq.gz,NanoMB1Rep1_R2.fastq.gz
```

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

### Outputs

| Directory         | Content                                                                                                                                                                               |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `reference/`      | The indexes that were built (`.fai`, `.dict`, `bwa/`).                                                                                                                                |
| `mapped_bams/`    | `<sample>_mapped_merged.bam`: the aligned BAM with the UMI tags.                                                                                                                      |
| `grouped_by_umi/` | The BAM grouped by UMI and the family size histogram.                                                                                                                                 |
| `consensus/`      | The duplex consensus reads: unmapped, mapped, filtered (`_filtered.bam`) and clipped for variant calling (`_filtered_clipped.bam`).                                                   |
| `variants/`       | `<sample>.vcf` from nvc and `<sample>_filtered.vcf`.                                                                                                                                  |
| `QC/`             | FastQC, FastQ Screen, MarkDuplicates metrics, samtools, Qualimap, `duplex_metrics/` and, with `--output_read_info`, `read_info/`. The consensus reads have their own `QC/consensus/`. |
| `QC/multiqc/`     | The MultiQC report (`multiqc_report.html`).                                                                                                                                           |
| `pipeline_info/`  | Tool versions (`software_versions.yml`), the execution report, timeline, trace and DAG.                                                                                               |

### Development

The tests use [nf-test](https://www.nf-test.com): `nf-test test --profile test,docker`. The FastQ Screen test needs a
database, so run it with `FASTQ_SCREEN_DB` set to the `FastQ_Screen_Genomes` folder made by
`tests/scripts/build_fastq_screen_db.sh` (bowtie2 is in `tests/scripts/environment.yml`). A GitHub Actions workflow
(`.github/workflows/nf-test.yml`) builds that database and runs all tests with docker on pushes and pull requests to `main`
and `dev`. The linters (prettier, editorconfig,
`nextflow lint`) run as [pre-commit](https://pre-commit.com) hooks. `plans/nextflow-rewrite.md` is the plan of the rewrite and
`CHANGELOG.md` lists the changes with respect to the Snakemake version.

## Snakemake version

The rest of this document describes the Snakemake pipeline (release 0.1.0).

## Prerequisites

Before running the pipeline you will need:

1. **Conda / Mamba** — all pipeline tools are installed automatically into
   isolated conda environments.
   [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) is
   recommended.
2. **Snakemake** — see [Installation](#installation) below.
3. **A BWA-indexed reference genome** — the reference FASTA specified in
   `config/config.yaml` (`ref`) must be indexed before running the pipeline:
   `bash bwa index ref/your_reference.fasta `
4. **Bowtie2 indexes for FastQ Screen** — the FastQ Screen config
   (`config/fastq_screen.conf`) must be configured to point to pre-built
   Bowtie2 indexes for any genomes you wish to screen against. See [FastQ
   Screen Configuration](#fastq-screen-configuration) below.

## Installation

Install Snakemake into a dedicated conda environment:

```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake=7.0.0-0
conda activate snakemake
```

Clone the repository:

```bash
git clone https://github.com/WEHIGenomicsRnD/duplex-seq-pipeline.git
cd duplex-seq-pipeline
```

## Input

Place paired-end FASTQ files in the `fastq/` directory. Files **must** follow the naming convention:

```
fastq/{sample}_R1.fastq.gz
fastq/{sample}_R2.fastq.gz
```

The pipeline automatically detects all samples from files matching this pattern
in the `fastq/` directory — no sample sheet is required.

## Configuration

Edit `config/config.yaml` before running. The key settings are:

| Parameter            | Description                                                                                                                                                                                                                                               |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ref`                | Path to the reference genome FASTA. Must be BWA-indexed (see [Prerequisites](#prerequisites)). Note that the directory containing the ref needs to be writeable to create the sequence dictionary file.                                                   |
| `read_structure`     | UMI read structure for each strand (e.g. `3M2S146T 3M2S146T`). See fgbio's [ExtractUmisFromBam](https://fulcrumgenomics.github.io/fgbio/tools/latest/ExtractUmisFromBam.html) for details.                                                                |
| `umis`               | Path to a text file of known UMI sequences (one per line). Leave empty to skip UMI correction.                                                                                                                                                            |
| `duplex_metrics`     | Metrics to compute with [calc-duplex-metrics](https://github.com/WEHIGenomicsRnD/calculate-duplex-metrics). Comma-separated list of individual metrics (`frac_singletons`, `efficiency`, `drop_out_rate`) and/or groups (`gc`, `family`). Default: `all`. |
| `variant_caller`     | Variant caller to use: `nvc` (naive variant caller) or `varscan`.                                                                                                                                                                                         |
| `min_reads`          | Minimum reads supporting a consensus base (fgbio `FilterConsensusReads --min-reads`).                                                                                                                                                                     |
| `FastqScreen_config` | Path to the FastQ Screen config file.                                                                                                                                                                                                                     |
| `cluster_config`     | Path to the SLURM cluster config YAML (per-rule resource settings). Default: `config/cluster.yaml`. Use `config/cluster_bigmem.yaml` for high-memory nodes.                                                                                               |
| `target_bed`         | For experiments with capture regions, you can provide a path in bed format to calculate the on-target rate for raw and duplex reads.                                                                                                                      |

All other parameters are documented inline in `config/config.yaml`.

### FastQ Screen Configuration

Edit `config/fastq_screen.conf` to add the genomes you want to screen against.
Each database entry must point to a pre-built **Bowtie2** index and be tagged
with `BOWTIE2`:

```
DATABASE    Human    /path/to/indexes/Human/GRCh38/Homo_sapiens.GRCh38  BOWTIE2
DATABASE    Mouse    /path/to/indexes/Mouse/GRCm38/Mus_musculus.GRCm38  BOWTIE2
DATABASE    Ecoli    /path/to/indexes/Ecoli/Ecoli                       BOWTIE2
```

Bowtie2 indexes can be built with:

```bash
bowtie2-build ref/your_reference.fasta /path/to/indexes/your_reference
```

## Running

You may be preferable to first run the pipeline to generate all the required conda environments:

```bash
snakemake --use-conda --conda-create-envs-only --cores 1
```

These may be a source of errors when running the pipeline, so a running this
first can help to make sure all envs are created successfully.

**Locally:**

```bash
conda activate snakemake
snakemake --use-conda --conda-frontend mamba --cores 4
```

**On a SLURM cluster:**

```bash
conda activate snakemake
snakemake --use-conda --conda-frontend mamba --profile slurm --jobs 8 --cores 24
```

## Outputs

All results are written under `results/`. The most relevant outputs are:

| Path                                                              | Description                                                                                                                                                                                        |
| ----------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `results/QC/multiQC/multiqc_report.html`                          | MultiQC report for pre-consensus reads, including FastQC and FastQ Screen metrics.                                                                                                                 |
| `results/QC/consensus/multiQC/multiqc_report.html`                | MultiQC report for post-consensus reads, including duplex-specific metrics from `results/QC/duplex_metrics/` via the [`duplex-multiqc`](https://github.com/WEHIGenomicsRnD/duplex-multiqc) plugin. |
| `results/QC/duplex_metrics/{sample}_metrics.csv`                  | Duplex metrics (efficiency, family statistics, GC bias, etc.) computed by [calc-duplex-metrics](https://github.com/WEHIGenomicsRnD/calculate-duplex-metrics).                                      |
| `results/consensus/{sample}_consensus_mapped_merged_filtered.bam` | Final duplex consensus read alignments (mapped, merged, and filtered).                                                                                                                             |
| `results/variants/{sample}_filtered.vcf`                          | Filtered variant calls.                                                                                                                                                                            |
