# Parity check against the Snakemake pipeline

`compare_parity.py` compares the outputs of the Snakemake pipeline (release 0.1.0) and of the Nextflow pipeline on the
same input, one line per artifact (BAMs read by read, consensus reads as multisets, UMI families as partitions, VCFs,
duplex metrics, samtools, Picard, Qualimap, FastQC, FastQ Screen). It needs `samtools` on the PATH.

```bash
python3 tests/parity/compare_parity.py <snakemake results dir> <nextflow outdir> [sample]
```

The results of the run on the test data are below (also in the notes of `CHANGELOG.md`).

## How the run was set up

1. Extract the Snakemake pipeline of the tag: `git archive 0.1.0 | tar -x -C snakemake_0.1.0`.
2. Give it the same input. Snakemake takes the sample name from the FASTQ file name before `_R1`, so link the test FASTQ
   files as `fastq/NanoMB1Rep1_R1.fastq.gz` and `_R2` (the Nextflow samplesheet calls the sample `NanoMB1Rep1`). Use the
   decompressed `tests/data` reference (`ref:` in `config/config.yaml`), build its bwa index with the bwa of the
   Snakemake environment (the pipeline expects an indexed reference), and set `target_bed` to `tests/data/target.bed`
   and `run_qualimap: true`.
3. FastQ Screen: the Snakemake `fastq_screen.conf` has no database, so point `FastqScreen_config` at a conf with one
   `DATABASE` line for the index made by `tests/scripts/build_fastq_screen_db.sh` (use an absolute path).
4. Run Snakemake with `--use-conda` on the same machine, then the Nextflow pipeline with
   `-profile test,<engine> --run_qualimap --fastq_screen_db <FastQ_Screen_Genomes>`, once with the default and once with
   `--markdup_duplex_umi false` (which marks duplicates per strand, as Snakemake does).

On a small machine the Snakemake rules had to be adapted (none of this changes what the tools compute): memory values
in a copy of `config/cluster.yaml` were capped, because Picard with `-Xmx90000m` was killed; the Java temp directory was
set with `JAVA_TOOL_OPTIONS=-Djava.io.tmpdir=...` (FastQC fails without a `/tmp`); and conda was pointed at a config with
only the conda-forge and bioconda channels. Snakemake installs `calculate-duplex-metrics` from the repository at run
time, so its version depends on the day of the run.

## Results on the test data

Both pipelines were run on the 1 kb test window (sample `NanoMB1Rep1`): Snakemake 8.11.3 with its own conda environments
and the Nextflow pipeline with apptainer, with `--markdup_duplex_umi false` and with the default. The script reports 21
comparisons, of which 5 are identical as they stand; the rest differ for the reasons below. The default run differs from
the per-strand run only in the Picard duplicate metrics and in the duplicate counts of samtools and Qualimap, which
follow from them.

**Identical**

- The aligned BAM, in all 86,593 records (position, MAPQ, CIGAR, mate fields, sequence, qualities and UMI tags), so bwa
  0.7.17 and 0.7.19 align identically. Only the `RG` tag differs (deliberate: the read group is the sample name, not
  `null`) and the `MQ` tag, which is missing on the 21 supplementary records in 0.1.0 (not isolated; probably GATK 4.7
  against Picard 2.27 MergeBamAlignment).
- The FastQ Screen results, the FastQC basic statistics, the (empty) filtered VCF and 19 of the 20 duplex metrics (to
  1e-9).
- The Picard MarkDuplicates metrics per strand (35,443 pair duplicates, 155 unpaired, 0.82291) apart from the optical
  duplicates.

**Different, with the cause.** Each cause was tested by running the tools of both versions on the same input.

| Artifact                  | Snakemake 0.1.0                         | Nextflow                       | Cause                                                                                                                                                                                                                                                                                                            |
| ------------------------- | --------------------------------------- | ------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| UMI grouping              | 85,572 records, 4,254 molecule families | 86,304 records, 4,359 families | fgbio 3 also groups the 366 templates (732 reads) with an unmapped read (241) or mates on different contigs (125); fgbio 1.4.0 drops them. All 4,254 families of 0.1.0 are unchanged.                                                                                                                            |
| Consensus reads           | 4,690                                   | 4,694                          | The new grouping adds 8; the consensus caller of fgbio 3 gives a different sequence for about 16% of the reads (3,947 of 4,690 are identical) and different qualities for more, with the same mean base quality (87.8 against 87.9). Each pipeline's consensus was reproduced exactly with the standalone fgbio. |
| Filtered consensus reads  | 3,380                                   | 3,580                          | `FilterConsensusReads` gives the same result in both versions on the same input; the difference comes from the consensus calls.                                                                                                                                                                                  |
| nvc sites                 | 801                                     | 761 (697 shared)               | Every site in both is an N-only call, so no variant is called and both filtered VCFs are empty. The sites that differ are consensus bases that one version calls as N.                                                                                                                                           |
| Duplex metrics            | 20 metrics                              | 19 metrics                     | Snakemake installs `calculate-duplex-metrics` from the repository at run time (commit `39fc6c9`), which has `unpaired_families`; the 0.2.0 container does not (0.2.1 is meant to).                                                                                                                               |
| samtools stats            | 101 outward oriented pairs              | 0                              | samtools 1.14 against 1.23 and later; the same on either BAM. The duplicate counts differ because the statistics read the marked BAM.                                                                                                                                                                            |
| Picard optical duplicates | 388                                     | 390                            | Picard 2.27.4 against 3.5.0 (2.27.4 gives 388 on the Nextflow BAM too).                                                                                                                                                                                                                                          |
| Qualimap                  | duplication rate estimated by Qualimap  | flagged duplicates             | The duplication lines (flagged against estimated), two extra lines for overlapping read pairs (the nf-core module runs Qualimap with `--collect-overlap-pairs`), and for the consensus BAM the consequences of the consensus differences.                                                                        |

With the default `--markdup_duplex_umi` Picard counts duplicates per duplex molecule (37,747 pair duplicates, 0.876148
instead of 35,443 and 0.82291), so the duplicate figures differ from 0.1.0 by design.

fgbio 3.1.2 is kept (as in the nf-core modules), so the consensus, the filtered BAM and the nvc calls do not match 0.1.0
exactly. Exact parity would need fgbio 1.4.0 in the local and nf-core fgbio modules.
