# Duplex Sequencing Pipeline

A [Snakemake](https://snakemake.readthedocs.io) pipeline for processing duplex
sequencing data (e.g.
[NanoSeq](https://www.nature.com/articles/s41586-021-03477-4)) to produce
duplex consensus reads, variant calls, and QC metrics. The workflow follows the
IDT recommended approach for [processing sequence data with unique molecular
identifiers](https://plone.bcgsc.ca/services/solseq/duplex-umi-documents/idt_analysisguideline_varcall-umis-dupseqadapters/).

## Prerequisites

Before running the pipeline you will need:

1. **Conda / Mamba** — all pipeline tools are installed automatically into
   isolated conda environments.
   [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) is
   recommended.
2. **Snakemake** — see [Installation](#installation) below.
3. **A BWA-indexed reference genome** — the reference FASTA specified in
   `config/config.yaml` (`ref`) must be indexed before running the pipeline:
   ```bash bwa index ref/your_reference.fasta ```
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

| Parameter        | Description                                                                                    |
|------------------|------------------------------------------------------------------------------------------------|
| `ref`            | Path to the reference genome FASTA. Must be BWA-indexed (see [Prerequisites](#prerequisites)). Note that the directory containing the ref needs to be writeable to create the sequence dictionary file. |
| `read_structure` | UMI read structure for each strand (e.g. `3M2S146T 3M2S146T`). See fgbio's [ExtractUmisFromBam](https://fulcrumgenomics.github.io/fgbio/tools/latest/ExtractUmisFromBam.html) for details. |
| `umis`           | Path to a text file of known UMI sequences (one per line). Leave empty to skip UMI correction. |
| `duplex_metrics` | Metrics to compute with [calc-duplex-metrics](https://github.com/WEHIGenomicsRnD/calculate-duplex-metrics). Comma-separated list of individual metrics (`frac_singletons`, `efficiency`, `drop_out_rate`) and/or groups (`gc`, `family`). Default: `all`. |
| `variant_caller` | Variant caller to use: `nvc` (naive variant caller) or `varscan`.                              |
| `min_reads`      | Minimum reads supporting a consensus base (fgbio `FilterConsensusReads --min-reads`).          |
| `FastqScreen_config` | Path to the FastQ Screen config file.                                                      |
| `cluster_config` | Path to the SLURM cluster config YAML (per-rule resource settings). Default: `config/cluster.yaml`. Use `config/cluster_bigmem.yaml` for high-memory nodes. |
| `target_bed`     | For experiments with capture regions, you can provide a path in bed format to calculate the on-target rate for raw and duplex reads. |

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

| Path                                               | Description                                                                        |
|----------------------------------------------------|------------------------------------------------------------------------------------|
| `results/QC/multiQC/multiqc_report.html`           | MultiQC report for pre-consensus reads, including FastQC and FastQ Screen metrics. |
| `results/QC/consensus/multiQC/multiqc_report.html` | MultiQC report for post-consensus reads, including duplex-specific metrics from `results/QC/duplex_metrics/` via the [`duplex-multiqc`](https://github.com/WEHIGenomicsRnD/duplex-multiqc) plugin. |
| `results/QC/duplex_metrics/{sample}_metrics.csv`   | Duplex metrics (efficiency, family statistics, GC bias, etc.) computed by [calc-duplex-metrics](https://github.com/WEHIGenomicsRnD/calculate-duplex-metrics). |
| `results/consensus/{sample}_consensus_mapped_merged_filtered.bam` | Final duplex consensus read alignments (mapped, merged, and filtered). |
| `results/variants/{sample}_filtered.vcf`           | Filtered variant calls.                                                            |
