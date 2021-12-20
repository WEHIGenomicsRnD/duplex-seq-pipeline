# Duplex Sequencing Pipeline

A [snakemake](https://snakemake.readthedocs.io) pipeline for obtaining duplex consensus reads generated using a duplex sequencing protocol (e.g. [NanoSeq](https://www.nature.com/articles/s41586-021-03477-4). The pipeline is similar to the recommend workflow from IDT on [processing sequence data with unique molecular
identifiers](https://plone.bcgsc.ca/services/solseq/duplex-umi-documents/idt_analysisguideline_varcall-umis-dupseqadapters/). QC is also performed, both pre- and post-duplex consensus calling. 

### Installation ###

The only prerequisite is [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). To install snakemake, you will need to install a Conda-based Python3 distribution. For this, [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) is recommended. Once mamba is installed, snakemake can be installed like so:

```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Now activate the snakemake environment (you'll have to do this every time you want to run the pipeline):

```
conda activate snakemake
```

Now clone the repository:

```
git clone https://github.com/WEHIGenomicsRnD/duplex-seq-pipeline.git
cd duplex-seq-pipeline
```

### Configuration ###

The configuration file is found under `config/config.yaml` and the config file for FastQ Screen is found under `config/fastq_screen.conf`. Please carefully go through these settings. The main settings to consider will be

- `read_structure` -- ensure that this matches the UMI design of your experiment. Refer to fgbio's [ExtractUmisFromBam](https://fulcrumgenomics.github.io/fgbio/tools/latest/ExtractUmisFromBam.html) for details on how to set this parameter. 
- `umis` -- if your UMIs are known (non-random), you can add a path to a UMIs text file (one UMI per line). Specifying a file path for this parameter will trigger a UMI correction step.
- `ref` -- ensure you have downloaded and specified the correct reference for your data.

### Running ###

Run the pipeline as follows:

```
conda activate snakemake
snakemake --use-conda --conda-frontend mamba --cores 1
```

If you want to submit your jobs to the cluster using SLURM, use the following to run the pipeline:

```
conda activate snakemake
snakemake --use-conda --conda-frontend mamba --profile slurm --jobs 8 --cores 24
```

The pipeline will generate all results under a `results` directory. The most relevant directories are:

- `results/QC/multiQC` -- contains the multiQC report for pre-consensus call reads. Also contains [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) metrics on the raw reads.
- `results/QC/consensus/multiQC` -- contains the multiQC report generated on reads after duplex consensus calling.
- ` results/consensus/{sample}__mapped_merged_filtered_clipped.bam` -- contains the mapped, filtered and clipped consensus reads. These should be your "final" read alignments for duplex consensus reads.
