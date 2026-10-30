# Plan: Nextflow rewrite of duplex-seq-pipeline (fastquorum-style, lightweight)

## Context

The Snakemake pipeline (`workflow/`) is being ported to Nextflow. Commit 087651c added a first attempt
in `nf-pipeline/`. It cannot run, and it is based on an older version of the Snakemake pipeline.

**Problems in the current skeleton**

- The 17 nf-core modules it includes are not installed.
- There is no `base.config`, schema, test data, versions output or `.fai` file.
- MultiQC is called with `Channel.empty()` inputs, so it never runs.
- The calls to the nf-core fgbio and bwa modules don't match their signatures, and the consensus parameters never reach the tools.
- `NVC` and `CALC_DUPLEX_METRICS` have placeholder `mulled-v2-XXX` containers.
- The `PREPROCESS` and `CONSENSUS` subworkflows are dead code.

**Snakemake changes the skeleton is missing**

| Commit  | Change                                                                                                                                                                                                                 |
| ------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 2960c23 | `run_qualimap` param (default false); qualimap becomes optional                                                                                                                                                        |
| c1550c2 | `get_read_info.py` rewritten: mate pairing, stricter filters, `epos` column, `sort -k 2,7` / `cut -f 2-7`                                                                                                              |
| e0a2a3f | read info merged into one step, with an empty-output check                                                                                                                                                             |
| 277c506 | `target_bed` param, passed as `--target_bed` to calc-duplex-metrics                                                                                                                                                    |
| d71e733 | `output_read_info` param (default false); calc-duplex-metrics reads the BAM by default (`--input_format bam --sort_mem --cores --bam_chunk_size`); GenomicAlignments added; duplex-multiqc plugin pinned at `41ab7fbc` |

**Goal**

A Nextflow pipeline that:

- is laid out like nf-core/fastquorum and follows nf-core best practice (meta maps, `ext.args`, versions, nf-schema validation, nf-test);
- leaves out nf-core-only extras (issue/PR templates, `.nf-core.yml`, ro-crate, igenomes, email, AWS/Tower CI, branch-protection and bot workflows);
- keeps the Snakemake defaults exactly;
- adds a samplesheet;
- produces one combined MultiQC report.

Containers are preferred over conda.

## Target layout (repo root = pipeline root; `nf-pipeline/` contents moved up — see Q4)

```
main.nf                         # entry: validate params, samplesheet, call workflow, write versions
workflows/duplexseq.nf          # main workflow
subworkflows/local/
  prepare_reference/            # samtools faidx, picard CreateSequenceDictionary, bwa index (skip if provided)
  fastq_align_umi/              # FastqToSam -> ExtractUmis -> [CorrectUmis] -> SamToFastq -> bwa mem -> MergeBamAlignment
  bam_align/                    # SamToFastq -> bwa mem -> MergeBamAlignment -> samtools index (shared by fastq_align_umi and consensus)
  consensus/                    # GroupReadsByUmi -> CallDuplex -> realign (same align steps) -> FilterConsensusReads
  qc/                           # fastqc, fastq_screen, markdup, samtools stats/flagstat, [qualimap], calc_duplex_metrics, [read_info]
  variant_calling/              # ClipBam -> nvc|varscan -> filter_vcf
  utils_pipeline/               # lightweight init: validateParameters, paramsSummaryLog, samplesheetToList (nf-schema only)
modules/nf-core/                # installed via `nf-core modules install` (tracked in modules.json)
modules/local/                  # fgbio extractumis/correctumis/clipbam, get_read_info, calc_duplex_metrics, nvc, varscan, filter_vcf, picard samtofastq (if no nf-core)
bin/                            # get_read_info.py, format_read_info.py, filter_vcf.py (synced from workflow/scripts at HEAD)
conf/base.config, conf/modules.config, conf/test.config
assets/schema_input.json, assets/multiqc_config.yml, assets/samplesheet.csv
nextflow.config, nextflow_schema.json, modules.json
tests/ (nf-test: default.nf.test + snapshot), nf-test.config, tests/.nftignore
.github/workflows/{lint.yml,nf-test.yml}
.editorconfig, .prettierrc.yml, .prettierignore, .pre-commit-config.yaml
README.md, docs/usage.md, docs/output.md, CHANGELOG.md
```

Removed from fastquorum's template: `.nf-core.yml`, ro-crate, `tower.yml`, `.devcontainer`, issue/PR templates,
`.dockstore.yml`, email templates, igenomes configs, methods-description template, `utils_nfcore_pipeline`,
and the `awstest`, `awsfulltest`, `branch`, `fix-linting`, `pr-comment`, `template-version-comment`,
`release-announcements`, `clean-up` and `download_pipeline` workflows.

## Analysis: Picard path vs fastquorum alignment path

**fastquorum path:** `fgbio FastqToBam` → `samtools fastq | bwa mem -K 150000000 -Y | fgbio ZipperBams | samtools sort --template-coordinate`

**Snakemake path:** `FastqToSam` → `ExtractUmisFromBam` → `SamToFastq` → `bwa mem -p` → `MergeBamAlignment`

Differences that change outputs:

1. **UMI tags (the most important one).**
   - `ExtractUmisFromBam --molecular-index-tags ZA ZB --single-tag RX` writes `ZA`/`ZB` as well as `RX`. `FastqToBam` writes only `RX` (`AAA-CCC`).
   - calc-duplex-metrics BAM mode (`R/read_info_from_bam.R`, 39fc6c9) and `get_read_info.py` both build the minus-strand UMI as `ZB-ZA`, and fetch `RX/ZA/ZB`.
   - With FastqToBam, every minus-strand pair would get UMI `NA`, so duplex pairing (x/y) and all duplex metrics would break unless we add a step to split `RX` into ZA/ZB or patch the R package.
2. **MergeBamAlignment defaults vs ZipperBams.**
   - MBA (Picard 2.27) defaults to `CLIP_ADAPTERS=true` and `CLIP_OVERLAPPING_READS=true`. It soft-clips read-through and overlap past the mate's start, and `INCLUDE_SECONDARY_ALIGNMENTS=true` keeps secondaries.
   - ZipperBams leaves bwa's alignments as they are.
   - Effects:
     - `epos` (CIGAR reference end) changes for short-insert pairs, which shifts read-bundle keys, so drop-out rate, efficiency and family counts move slightly.
     - samtools stats "bases mapped (cigar)" and error rate shift.
     - GroupReadsByUmi uses unclipped 5′ positions, so grouping is barely affected.
3. **bwa flags.** `-K 150000000` makes results independent of thread count, which is better for reproducibility. `-Y` soft-clips supplementaries; we filter supplementaries out anyway. Minor.
4. **Sort order.** fastquorum produces template-coordinate output. MarkDuplicates, samtools stats, qualimap and calc-duplex-metrics need coordinate order (calc-duplex-metrics pairs mates in chunks and assumes coordinate order), so we'd need an extra sort and BAM. Little net saving.
5. **Tool versions, which matter whichever path we choose.** Snakemake pins **fgbio 1.4.0, picard 2.27.4, bwa 0.7.17, samtools 1.14, fastqc 0.11.8, fastq-screen 0.14.0, qualimap 2.2.2, varscan 2.4.0, nvc 0.0.4**. nf-core modules use fgbio 2.5.x and picard 3.x.
   - Exact parity needs local modules (or nf-core modules with container overrides) pinned to those versions. Biocontainers exist for all of them.
   - fgbio 2.x changes GroupReadsByUmi/consensus behaviour in places, for example `--allow-inter-contig` and the sort-order requirements.

**Decisions (user-confirmed):**

- **Alignment:** keep the Picard path (FastqToSam/ExtractUmis/MBA). `-K 150000000` can be added to bwa through `ext.args` but stays off by default.
- **Versions:** upgrade to current nf-core module versions (fgbio 3.1.2, picard 3.5.0, GATK 4.7.0.0, bwa 0.7.19, samtools 1.22–1.25; these are what the installed nf-core modules use). Validate against Snakemake on the test data and record the differences in CHANGELOG/docs.
  - Check after the upgrade: fgbio 2.x GroupReadsByUmi still accepts a coordinate-sorted MBA BAM (it sorts internally if needed).
  - Check after the upgrade: `ExtractUmisFromBam` still exists and keeps `--molecular-index-tags`.

## Containers needed

| Tool                                                                             | Container source                                                                                                                                                                                                                                                                                 |
| -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| fastqc, fastq-screen (+bowtie2), picard, fgbio, bwa+samtools, samtools, qualimap | nf-core modules (Seqera/biocontainers); pinned biocontainer tags if we keep the Snakemake versions                                                                                                                                                                                               |
| varscan + samtools, nvc + samtools                                               | multi-tool: Seqera Containers built from the bioconda spec (nvc 0.0.4 is on bioconda)                                                                                                                                                                                                            |
| python/pysam/pandas (`get_read_info.py`, `format_read_info.py`, `filter_vcf.py`) | Seqera Containers built from the conda spec                                                                                                                                                                                                                                                      |
| calc-duplex-metrics                                                              | **exists**: `ghcr.io/wehigenomicsrnd/calculate-duplex-metrics:<release tag>` (build-docker.yml in that repo). Built on release. Latest is `0.2.0`; check that it includes BAM mode and 39fc6c9, otherwise cut `0.2.1`. `calc-duplex-metrics` is on PATH at `/usr/local/bin`. No runtime install. |
| MultiQC + duplex-multiqc plugin                                                  | **needs building**: add a Dockerfile (`FROM` the Seqera multiqc 1.3x image, `pip install duplex-multiqc@<tag>`) and a GHCR build workflow in the duplex-multiqc repo. Conda profile: `multiqc` + pip install from git.                                                                           |

## Steps

### 1. Scaffold and config

- Move `nf-pipeline/` to the repo root. The Snakemake code stays under `workflow/` until parity is confirmed.
- **`nextflow.config`:**
  - `manifest` (`nextflowVersion '!>=25.04.0'`) and plugin `nf-schema@2.5.1`;
  - `docker.registry` and `singularity.registry = 'quay.io'`;
  - `process.resourceLimits` replaces the unused `max_*` params;
  - profiles: `docker`, `singularity`, `apptainer`, `conda`, `mamba`, `wave`, `test`, `debug`;
  - `includeConfig` for `conf/base.config` and `conf/modules.config`;
  - `overwrite = true` on the timeline, report, trace and dag outputs.
- **`conf/base.config`:** copy fastquorum's labels (`process_single` … `process_high_memory`) and its retry strategy.
  - Picard, fgbio and calc_duplex_metrics need large memory, so map the Snakemake `cluster.yaml` defaults onto `withName` overrides in `base.config`:
    - picard: 96 GB
    - bwa: 16 cpus, 36 GB
    - calc_duplex_metrics: 16 cpus, 128 GB
    - fgbio: 64 GB
  - Java heap is `-Xmx` at about 90% of `task.memory`, matching how `java_mem` was used.

### 2. Params, identical to `config/config.yaml`

- Same names and values as Snakemake:
  - `read_structure='3M2S146T 3M2S146T'`, `molecular_index_tags='ZA ZB'`, `single_tag='RX'`
  - `umis=''`, `max_mismatches=3`, `min_distance=1`
  - `min_mapq=20`, `min_input_baseq=30`, `error_rate_pre_umi=45`, `error_rate_post_umi=40`
  - `min_reads='4 2 2'`, `max_read_error_rate=0.05`, `max_base_error_rate=0.1`, `min_base_quality=50`, `max_no_call_fraction=0.05`
  - `duplex_metrics='all'`, `target_bed=null`, `output_read_info=false`, `run_qualimap=false`
  - `variant_caller='nvc'`, `min_coverage=4`, `min_reads2=2`, `min_vaf=0.01`, `ploidy=1`
- Resource-level values also become params:
  - `sort_mem='4G'`, `bam_chunk_size=2000000`
  - `read_info_buffer_size='110G'`
  - `markdup_max_records_in_ram=100000`
- Renamed:
  - `ref` → `fasta`, with optional `fasta_fai`, `dict` and `bwa` for prebuilt indexes;
  - `FastqScreen_config` → `fastq_screen_config` (FastQ Screen is skipped when unset).
- Kept hard-coded, as in Snakemake, but set as `ext.args` in `modules.config` so they can be overridden:
  - GroupReadsByUmi `--strategy=paired --edits=1`
  - MergeBamAlignment `--SORT_ORDER coordinate --ALIGNER_PROPER_PAIR_FLAGS true --MAX_GAPS -1 --ORIENTATIONS FR` (the index comes from `samtools index`)
  - ClipBam `--clipping-mode=Hard --clip-overlapping-reads=true`
  - MarkDuplicates `BARCODE_TAG=RX`
  - nvc `--min_base_quality=0`
  - varscan `--min-avg-qual 0 --strand-filter 0 --p-value 1`
- `nextflow_schema.json` is generated with `nf-core pipelines schema build` as a one-off local tool; nf-core is not used in CI. Params are grouped by step, mirroring the fastquorum schema.

### 3. Samplesheet (`assets/schema_input.json`)

- Columns: `sample, fastq_1, fastq_2, [read_structure], [umi_file]`. The two optional columns override `params.read_structure` and `params.umis` per sample.
- Parsed with `samplesheetToList`. Multi-lane merging is out of scope for v1; the Snakemake pipeline doesn't have it either.

### 4. Modules (containers first)

**nf-core modules** (Seqera/biocontainers images):

| Area                   | Modules                                                                                                                             |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| QC                     | `fastqc`, `fastqscreen/fastqscreen`, `qualimap/bamqc`, `multiqc`                                                                    |
| Picard                 | `picard/createsequencedictionary`, `picard/fastqtosam`, `gatk4/mergebamalignment`, `gatk4/samtofastq`, `picard/markduplicates`      |
| Alignment and samtools | `bwa/index`, `bwa/mem` (sort=false, two FASTQs, no `-p`), `samtools/faidx`, `samtools/stats`, `samtools/flagstat`, `samtools/index` |
| fgbio                  | `fgbio/groupreadsbyumi`, `fgbio/callduplexconsensusreads`, `fgbio/filterconsensusreads`                                             |

- **Wiring:** check each signature after install and pass params through `val` inputs or `ext.args`, so that every Snakemake argument reaches its tool (this fixes the skeleton's mismatched calls).
- **Local modules** (meta map, `versions.yml` or topic `versions`, `stub:` block, `ext.args`/`ext.prefix`, no `params.*` inside the module):
  - **`FGBIO_EXTRACTUMISFROMBAM`, `FGBIO_CORRECTUMIS`, `FGBIO_CLIPBAM`:** Seqera fgbio container (or an nf-core module if one exists).
  - **SamToFastq / MergeBamAlignment:** not in `picard/` but in nf-core as `gatk4/samtofastq` and `gatk4/mergebamalignment` (GATK 4.7 bundles Picard). No local `PICARD_SAMTOFASTQ` is needed. `gatk4/samtofastq` writes `_1`/`_2` FASTQs (it cannot interleave), so `bwa mem` reads two files without `-p`; the alignments are the same. MergeBamAlignment arguments use the new `--NAME value` Picard syntax, and `CREATE_INDEX` is replaced by a `samtools/index` step.
  - **`GET_READ_INFO`:** one process, ported from the current `qc.smk` rule:
    - runs `get_read_info.py | sort -k 2,7 -k 1,1 -u | cut -f 2-7 | uniq --count`, then `format_read_info.py`;
    - fails when the output is empty;
    - `bin/` scripts replaced with the HEAD versions from `workflow/scripts/`;
    - container: Seqera image built from the conda spec (python, pysam, pandas).
  - **`CALC_DUPLEX_METRICS`:**
    - BAM mode by default, read-info mode when `output_read_info` is set;
    - optional `--target_bed`;
    - Snakemake's `--metrics`, `--sort_mem`, `--cores ${task.cpus}` and `--bam_chunk_size` options;
    - no runtime `install_github`; the container comes from Q2.
  - **`NVC`, `VARSCAN`, `FILTER_VCF`:**
    - containers are Seqera images built from the bioconda specs (nvc, varscan + samtools, python);
    - placeholder `mulled-v2-XXX` images removed;
    - the varscan branch also gets `FILTER_VCF` applied or skipped consistently (fixes the Snakemake bug where `rule all` can't be satisfied).
- **Picard input order:**
  - feed MergeBamAlignment a queryname/unsorted aligned BAM, not a coordinate-sorted one;
  - set `ext.prefix` so the unmapped and aligned BAMs don't share a filename in the consensus merge;
  - emit the `.bai`.

### 5. Single combined MultiQC report

- One `MULTIQC` call. Inputs:
  - FastQC zips, FastQ Screen txt, MarkDuplicates metrics, raw samtools stats/flagstat, optional qualimap;
  - consensus samtools stats/flagstat and optional consensus qualimap;
  - `*_metrics.csv` from calc_duplex_metrics;
  - software versions YAML and the params summary.
- **Naming:** consensus outputs get `ext.prefix = "${meta.id}_consensus"`; raw outputs keep `${meta.id}` (FastQC gives `${meta.id}_R1`/`_R2`).
- **`assets/multiqc_config.yml`:**
  - `table_sample_merge: { "R1": "_R1", "R2": "_R2", "Consensus": "_consensus" }` groups R1, R2, raw and consensus as sub-rows under one sample in the General Statistics table (needs MultiQC ≥ 1.22);
  - `extra_fn_clean_exts` is limited to `_mapped_merged`, so `_consensus` is not stripped and raw and consensus samples don't collide in the plots;
  - `report_section_order` puts duplex metrics first;
  - `software_versions` section.
- **Plugin:** pip `duplex-multiqc` on a current MultiQC (1.3x). The upstream plugin needs updating for the new MultiQC API (user will do this in the plugin repo). The container comes from Q2.

### 6. Test data

**Supplied (in `tests/data/`):**

| File                                        | Content              |
| ------------------------------------------- | -------------------- |
| `NanoMB1Rep1_HJK2GDSX3_sampled_R1.fastq.gz` | 43,286 reads, 3.5 MB |
| `NanoMB1Rep1_HJK2GDSX3_sampled_R2.fastq.gz` | 43,286 reads, 3.6 MB |

- **What it is:** a 1 kb window of the E. coli ATCC 10798 duplex sequencing data (43,286 read pairs overlapping a 1 kb window, extracted by `extract_region_reads.sh` in the Snakemake repo; all reads map to contig `NARG01000009.1` of the test reference), with every read pair in the window kept so UMI families stay at full depth.
  - This is the "take a region, not a random fraction" approach: random thinning (`samtools view -s`) would leave almost no duplexes with `min_reads '4 2 2'`.
- **Read structure:** NanoSeq, `3M2S146T 3M2S146T`.
- **Checked in the files:**
  - R1 and R2 have identical read names in the same order, with no duplicate names.
  - Every read is **150 bp**.
  - The first 3 bases (the UMI) are diverse (81 distinct values including N-containing ones), so UMIs are intact and not yet extracted.
  - Bases 4–5 (the skipped `2S`) are strongly non-random (`GT`, `CT`, `TC`, `TA` dominate), consistent with the NanoSeq spacer.
- **Not yet checked (needs bwa/samtools, see Readiness):** that the window really is about 1 kb, and how many pairs end up as duplex families.
- **Expected scale:** about 43k pairs over 1 kb is roughly 40x depth per strand-pair position. The UMI space is only 64 × 64 = 4,096 combinations, so UMI collisions within the window are likely and will affect family counts. This is realistic for NanoSeq, but compare against Snakemake on the same data, not against theory.

**Still needed (not in the repo):**

1. **Reference FASTA** (done): `tests/data/GCA_002087985.1_ASM208798v1_genomic.fna.gz` (GenBank assembly ASM208798v1 of E. coli ATCC 10798, 64 contigs; **plain gzip, not bgzip**, so the workflow decompresses it with `GUNZIP_FASTA` before indexing). **The strain is ATCC 10798, not BL21; the plan's earlier BL21 references were wrong.** An earlier private copy of the reference (2 circular contigs) was dropped from the repo because its licence could not be confirmed; its contig names differ, so window coordinates from `extract_region_reads.sh` do not apply to the new reference.
   - Copy it (gzipped, about 1.4 MB) into `tests/data/` for the test profile; bwa indexing takes seconds, so no prebuilt index is committed.
   - The FASTQs were produced by aligning the full-genome reads to a reference and pulling pairs overlapping a 1 kb window. Mates of those pairs may map outside the window (insert sizes are mostly 237–353 bp, but some are much larger), so keep the whole genome rather than a slice. Against the new reference the read counts are identical (86,572 primary records, 99.72% mapped), with essentially all reads on `NARG01000009.1` (around 108,900–111,500).
2. **`target.bed`**: a sub-window inside the 1 kb on `NARG01000009.1` (exercises `--target_bed`). The exact window must be re-derived from read coverage on the new reference.
3. **`umis.txt`** (optional): the NanoSeq 3-mer whitelist, to exercise CorrectUmis. If there is no real whitelist, generate all 64 ACGT 3-mers; that is not a realistic test, because every UMI is then valid, but it does run the step. Without it, `--umis` is untested.
4. **Second sample** (optional): a different window, or the same FASTQs under a second sample name, to test multi-sample MultiQC merging. Same-FASTQ duplication is cheap and fine for a smoke test.

**Read length (resolved empirically).** `3M2S146T` consumes 151 bases but the test reads are 150 bp. I expected `ExtractUmisFromBam` to fail; it does not.

- Test: `samtools import` of the two FASTQs into an unmapped BAM, then `ExtractUmisFromBam -r <rs> -t ZA ZB -s RX` with **fgbio 1.4.0 and 2.5.21**, each with `3M2S146T 3M2S146T` and `3M2S145T 3M2S145T`.
- All four runs succeed and produce **byte-identical BAMs** (same md5): 86,572 records, templates trimmed to 145 bp, with `ZA`, `ZB` and `RX` set (for example `ZA:GTG ZB:ATA RX:GTG-ATA`). fgbio simply truncates the final fixed-length template segment when the read is one base short.
- **Decision:** keep the default `read_structure = '3M2S146T 3M2S146T'` in the test profile. No override and no re-extraction are needed. Real 151 bp data and this 150 bp data behave the same for UMI extraction.
- This also confirms for the version upgrade that `ExtractUmisFromBam` exists in fgbio 2.5.21 and keeps `--molecular-index-tags` (`-t`) and `--single-tag` (`-s`).

**Hosting:** the data lives in this repo under `tests/data/`, so `conf/test.config` uses relative paths (`${projectDir}/tests/data/...`) and no raw URLs. Total size is under 10 MB. `-profile test,docker` should finish in under 5 minutes.

**Test samplesheet:** `tests/data/samplesheet.csv` (committed, relative paths), columns `sample,fastq_1,fastq_2`. Sample name: `NanoMB1Rep1`.

### 7. Tests and CI (lightweight)

- **nf-test:** `tests/default.nf.test` snapshots the versions, published file names and stable files (VCF body, duplex metrics CSV), with BAMs and MultiQC HTML ignored in `.nftignore`. Plus a stub run test.
- **`.github/workflows/lint.yml`** runs pre-commit via prek with these hooks:
  - `prettier`
  - `editorconfig-checker` (`editorconfig-checker/editorconfig-checker.python`)
  - `trailing-whitespace` and `end-of-file-fixer`, excluding `modules/nf-core` and `subworkflows/nf-core`
  - `nextflow-lint` (`seqeralabs/nf-lint-pre-commit`), a syntax and strict-syntax check
  - `check-jsonschema` on `nextflow_schema.json` and `assets/schema_input.json`, validated against the JSON Schema 2020-12 meta-schema
  - a small `nextflow config -flat` / `nextflow run . --help` smoke step, which checks that the schema loads and every param is in the schema (nf-schema's `validateParameters` fails on unknown or undeclared params)
- **`.github/workflows/nf-test.yml`:** nf-test with `-profile test,docker` on `ubuntu-latest`, Nextflow pinned and `latest-everything`, runs on PRs.
- Copied from fastquorum and trimmed: `.editorconfig` (new), `.prettierrc.yml`, `.prettierignore`.

### 8. Docs

- README: quick start, samplesheet, params table mapped from the old config keys.
- `docs/output.md`.
- Note in README that Snakemake is deprecated, once parity is verified.

## Verification

1. **Static checks:** `prek run --all-files` (prettier, editorconfig, nextflow lint, schema check) and `nextflow run . --help`.
2. **Stub run:** `nextflow run . -profile test,docker -stub` to check that channels are wired correctly.
3. **Real run:** `nextflow run . -profile test,docker` on the 1 kb E. coli window in `tests/data/`. Every step should run and a single `multiqc/multiqc_report.html` should be produced, with R1/R2/raw/consensus grouped per sample in General Statistics and the duplex metrics section present.
4. **Parity with Snakemake:** run the Snakemake pipeline on the same test data and compare:
   - consensus BAM read counts (`samtools flagstat`)
   - `*_metrics.csv` (should be identical)
   - filtered VCF records
   - read-info output with `--output_read_info`
5. **Optional branches:** `--run_qualimap`, `--output_read_info`, `--umis`, `--target_bed`, `--variant_caller varscan`.
6. `nf-test test --profile docker` passes locally and in CI.

## Readiness check (before implementation)

Checked against this repo and the login node the session runs on.

**Have**

- The plan, with all design decisions made (Picard path, upgraded tool versions, containers first, one MultiQC).
- Snakemake source for every step (`workflow/rules`, `workflow/scripts`, `config/`), and the old skeleton in `nf-pipeline/` to mine for module logic.
- Paired test FASTQs (`tests/data/`), consistent with each other.

**Missing: needs you**

| #   | Item                                                               | Blocks                               | Notes                                                                                                                                                |
| --- | ------------------------------------------------------------------ | ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | ~~Reference FASTA and window coordinates~~                         |                                      | **Done**: `GCA_002087985.1` (ATCC 10798) FASTA in `tests/data/`; reads map to `NARG01000009.1`. Window coordinates on this reference still to derive |
| 2   | ~~Read length: 150 bp data vs `3M2S146T`~~                         |                                      | **Resolved**: fgbio truncates the last segment; no change needed (section 6)                                                                         |
| 3   | `umis.txt` whitelist (optional)                                    | CorrectUmis test only                | Can be skipped for session 1                                                                                                                         |
| 4   | calc-duplex-metrics image: does `0.2.0` have BAM mode and 39fc6c9? | `CALC_DUPLEX_METRICS` container tag  | Needs a check of that repo, or cut `0.2.1`                                                                                                           |
| 5   | duplex-multiqc container and the plugin update for MultiQC 1.3x    | Real MultiQC run with duplex metrics | You said you will do the plugin update; until then session 1 can use plain MultiQC and skip the duplex section                                       |

**Environment (re-checked; now usable, with workarounds)**

- Git worktree works. Modules work: `source /etc/profile.d/modules.sh; module load nextflow/25.10.0 openjdk/21.0.2 samtools/1.23.1 bwa/0.7.17 nf-core/2.14.1` (Nextflow 25.10.0, Java 21; `nf-core` 2.14.1 and `pre-commit` come with the nf-core module; `nf-test` and `prek` are not installed).
- **Network works, but needs a CA bundle** (`/etc/pki/tls/certs/ca-bundle.crt` is missing in the sandbox). Set `CURL_CA_BUNDLE`, `SSL_CERT_FILE` and `REQUESTS_CA_BUNDLE` to `/stornext/System/data/software/rhel/9/base/tools/miniconda3/latest/lib/python3.12/site-packages/certifi/cacert.pem`. Set each variable to the literal path, not to `$CURL_CA_BUNDLE` in the same `export` (it expands empty). Java downloads (Nextflow plugins) are untested.
- **Containers now run** (checked with fastqc, fgbio, calc-duplex-metrics and a full pipeline run). The sandbox needs these workarounds, which are kept out of the pipeline in the git-ignored `conf/sandbox.config` (`-c conf/sandbox.config -profile test,sbx_apptainer`):
  - the module-provided `~/.nextflow/config` forces `executor = 'slurm'` and singularity, which fail here; the sandbox config sets `process.executor = 'local'` and disables singularity;
  - `/var/apptainer/mnt/session` must exist; apptainer also needs a `root` entry in `/etc/passwd`;
  - `--no-mount bind-paths,sys,tmp` (as `apptainer.runOptions`; `APPTAINER_NO_MOUNT` does not reach the task because Nextflow starts apptainer with `env -`), plus `-B <dir>:/tmp` for tools that need a temp dir (the R-based calc-duplex-metrics does);
  - `apptainer.pullTimeout` must be set explicitly (`'30 min'`), otherwise Nextflow fails with a null `pullTimeout`;
  - pulls fail with an xattr error when `APPTAINER_TMPDIR` is on NFS; point `APPTAINER_TMPDIR` and `APPTAINER_CACHEDIR` at the sandbox tmpfs;
  - load only `nextflow/25.10.0` (`module purge` first): loading `openjdk` as well made the launcher fall back to Nextflow 23.10.0, which ignores user plugins.
  - Helper used this session (not in the repo): `/vast/scratch/users/cmero.ma/tmp/nfenv.sh` sets up the modules, CA bundle and apptainer environment.
  - Shared cache `/vast/projects/genomics_cache/singularity_cache/` (read-only) has fastqc 0.12.1, qualimap 2.3, samtools 1.23.1, multiqc 1.33 and `calculate-duplex-metrics-0.2.0` (which has BAM mode: `--input_format {rinfo,fgbio,bam}`, `--bam_chunk_size`, `--target_bed`); personal cache `/vast/scratch/users/cmero.ma/my_singularity_cache/` now also holds the Seqera images pulled by the pipeline.
- `nf-core` 2.14.1 (module) installs modules from nf-core/modules master, but it needs a pipeline marker: use a temporary `.nf-core.yml` (`repository_type: pipeline`, git-ignored) and an existing `modules/` directory.

**Decisions that affect session 1**

- nf-core modules can't be installed offline. Options: (a) you run `nf-core modules install ...` (list below) and commit `modules.json` and `modules/nf-core/`; (b) I write the nf-core modules by hand as local copies. I recommend (a): hand-copied nf-core modules drift, and `modules.json` tracking is part of the design.
- Fix the git worktree link, or tell me where the repo should be committed from.

## Session 1 plan

**Goal:** a complete, lint-clean Nextflow skeleton at the repo root, with the preprocessing half of the workflow (reference prep, FASTQ to merged aligned BAM) wired and, if the environment allows, run on the test data. Later sessions add consensus, QC/MultiQC and variant calling.

**Why this slice:** it exercises the riskiest early assumptions (Picard/fgbio path, `ZA`/`ZB` UMI tags, MBA input order, tool-version upgrade) and everything downstream depends on it. Nothing in it needs the missing duplex-multiqc container or the calc-duplex-metrics image.

**Part A: no tools needed (always done in session 1)**

1. Move `nf-pipeline/` to the repo root (Snakemake stays under `workflow/`); delete the dead `PREPROCESS`/`CONSENSUS` subworkflows and placeholder modules, keeping the logic worth reusing.
2. `nextflow.config`, `conf/base.config`, `conf/modules.config` (with every Snakemake default as `ext.args`), `conf/test.config`; params exactly as in section 2.
3. `assets/schema_input.json`, `tests/data/samplesheet.csv`, `nextflow_schema.json` written by hand (nf-core schema tooling isn't available) and checked against the JSON Schema meta-schema.
4. `subworkflows/local/utils_pipeline` (nf-schema validation, samplesheet parsing).
5. Local modules for the preprocessing half: `FGBIO_EXTRACTUMISFROMBAM`, `FGBIO_CORRECTUMIS`, `PICARD_SAMTOFASTQ` (if no nf-core one), each with `ext.args`/`ext.prefix`, `stub:` and versions.
6. Subworkflows `prepare_reference` and `fastq_align_umi`, and `workflows/duplexseq.nf` + `main.nf` calling them.
7. `bin/` scripts synced from `workflow/scripts` at HEAD.
8. `.editorconfig`, `.prettierrc.yml`, `.prettierignore`, `.pre-commit-config.yaml`, `lint.yml`.
9. README stub: quick start, samplesheet, params table.

**Part B: needs nextflow + network (done if available, otherwise handed over as a checklist)**

1. You or I run `nf-core modules install` for: `fastqc`, `fastqscreen/fastqscreen`, `picard/createsequencedictionary`, `picard/fastqtosam`, `picard/mergebamalignment`, `bwa/index`, `bwa/mem`, `samtools/faidx`. (The remaining modules are installed in session 2.)
2. Reconcile each installed module's inputs with the wiring in Part A; check the points flagged in the Decisions section (`ExtractUmisFromBam` and `--molecular-index-tags` still exist in fgbio 2.x; `bwa mem -p`; `MergeBamAlignment` arguments).
3. `nextflow run . -profile test,docker -stub`, then a real run to the merged BAM. Check that `ZA`/`ZB`/`RX` tags are present in the output BAM (`samtools view | head`) and that reads aligned to the reference.
4. (The 150 bp vs `3M2S146T` question is already resolved; see section 6.)

**Out of scope for session 1:** consensus calling, QC subworkflow, MultiQC, variant calling, nf-test snapshots, the full docs, parity comparison with Snakemake, containers for NVC/varscan/python.

**Done when**

- Part A files exist and `nextflow config` plus `nextflow run . --help` load cleanly (if Nextflow is available), and the schema files validate.
- If the environment allows Part B: a stub run passes and the real run produces an aligned BAM with UMI tags from the test data.
- Anything not run is listed explicitly in the end-of-session notes, not assumed to work.

**Session 2 and 3 (outline)**

- Session 2: `consensus` subworkflow (GroupReadsByUmi → CallDuplex → realign → FilterConsensusReads), the remaining nf-core installs, `qc` subworkflow without duplex metrics, single MultiQC.
- Session 3: `CALC_DUPLEX_METRICS`, `GET_READ_INFO`, `variant_calling`, containers for local tools, nf-test snapshots, CI, docs, Snakemake parity run.

## Session 1 results

**Done (Part A and Part B both completed)**

- Repo restructured: `nf-pipeline/` moved to the repo root (Snakemake stays under `workflow/`); dead subworkflows, placeholder modules and the old `main.nf`/`nextflow.config` removed (still in git history at 087651c if the logic is wanted for later modules).
- Config: `nextflow.config` (nf-schema 2.5.1, profiles, all params of section 2), `conf/base.config`, `conf/modules.config`, `conf/test.config`; schemas `nextflow_schema.json` (hand-written, cross-checked against the config: every param present, defaults equal) and `assets/schema_input.json`.
- Code: `main.nf`, `workflows/duplexseq.nf`, subworkflows `utils_pipeline`, `prepare_reference`, `bam_align`, `fastq_align_umi`; local modules `FGBIO_EXTRACTUMISFROMBAM` and `FGBIO_CORRECTUMIS` (fgbio 3.1.2, same container as the nf-core fgbio modules).
- nf-core modules installed and tracked in `modules.json`: `fastqc`, `fastqscreen/fastqscreen`, `picard/createsequencedictionary`, `picard/fastqtosam`, `gatk4/mergebamalignment`, `gatk4/samtofastq`, `bwa/index`, `bwa/mem`, `samtools/faidx`, `samtools/index`.
- Test data: reference FASTA and `tests/data/samplesheet.csv` added next to the FASTQs.
- Tooling: `.editorconfig`, `.prettierrc.yml`, `.prettierignore`, `.pre-commit-config.yaml`, `.github/workflows/lint.yml`, `.gitignore`; README gains a Nextflow section.

**Verified**

- `nextflow run . --help` works (plugin downloads); parameter validation rejects bad values (for example `--variant_caller foo`, a missing `--dict`).
- `nextflow lint` and `pre-commit run --all-files` (prettier, editorconfig, whitespace, nextflow-lint, metaschema check) pass.
- Stub run (`-profile test -stub`) passes through all processes, also with `--umis`, `--bwa`, `--dict` and `--fasta_fai` supplied (reference steps skipped).
- **Real run** with apptainer on the test data (`-profile test,sbx_apptainer`) completes: `mapped_bams/NanoMB1Rep1_mapped_merged.bam` (+ `.bai`) has 86,572 primary reads (all pairs), 99.7% mapped, coordinate-sorted, read group `null` / sample `NanoMB1Rep1`, templates 145 bp, and `ZA`/`ZB`/`RX` on every record (for example `ZA:TCA ZB:TAT RX:TCA-TAT`). `pipeline_info/software_versions.yml` lists the tool versions from the `versions` topic.
- Real run with `--umis` (all 64 3-mers, written to scratch only) runs `FGBIO_CORRECTUMIS` with fgbio 3.1.2. It drops 118 pairs, consistent with UMIs that contain `N` and cannot match any whitelist entry.

**Deviations from the plan (all recorded above)**

- `picard/mergebamalignment` does not exist; `gatk4/mergebamalignment` and `gatk4/samtofastq` are used (GATK 4.7.0.0). `samtofastq` cannot interleave, so `bwa mem` takes two FASTQs and no `-p`.
- MergeBamAlignment no longer sets `CREATE_INDEX`; `samtools/index` produces the `.bai`.
- fgbio is 3.1.2 (not 2.5.x): that is what nf-core now ships. `ExtractUmisFromBam` and `CorrectUmis` still work with the Snakemake arguments. fgbio 1.4.0 and 2.5.21 gave identical `ExtractUmisFromBam` output on the test data (section 6), and 3.1.2 behaved the same in the real run.
- Added `bam_align` as a shared subworkflow, so the consensus step can reuse the realignment.
- Samplesheet FASTQ and UMI-file paths are plain strings in `schema_input.json` (no `format: file-path`/`exists`), and `utils_pipeline` resolves relative paths against the samplesheet's directory: nf-schema resolves relative paths against the launch directory, so the committed test samplesheet would not have worked otherwise. Missing files are reported by `file(..., checkIfExists: true)`.
- Test config uses `${projectDir}/tests/data/...` (no URLs), as decided in section 6.
- `ext.args` for FastqToSam sets `--SAMPLE_NAME ${meta.id} --READ_GROUP_NAME null` (the module otherwise takes the sample name from `ext.prefix`).

**Open points and things to know for session 2**

- **Nextflow syntax:** `nextflow lint -format` rewrites `String path` to `path: String`, which Nextflow 25.10's default parser rejects (it needs the v2 parser, the default from 26.04). Don't use typed parameters or declarations while `nextflowVersion` is `!>=25.04.0`.
- **Not yet done:** no `meta.yml` for the local modules; no nf-test files; `nf-test.yml` workflow; docs beyond the README section; the per-module `publishDir` layout is only set up for the session 1 outputs.
- **Not yet compared with Snakemake:** the merged BAM has not been compared with `results/mapped_bams/*_mapped_merged.bam` from the Snakemake pipeline (tool versions differ: bwa 0.7.19 vs 0.7.17, picard 3.5/GATK 4.7 vs 2.27.4, samtools 1.22 vs 1.14). The merged BAM contains 21 supplementary records (bwa chimeric alignments); Snakemake keeps them too, and the QC steps filter on `-F 0x900`.
- **Naming:** `manifest.name` and the schema `$id`s use `mcmero/duplex-seq-pipeline` (the git `origin`), while the README clones `WEHIGenomicsRnD/duplex-seq-pipeline`. Decide which is canonical before the first release. Contributors are not listed in the manifest.
- **Git:** nothing is committed. The new files are marked intent-to-add (`git add -N`) so pre-commit could see them, and the index still holds the earlier `nf-pipeline/` renames.
- **Resources:** `conf/base.config` uses the Snakemake `cluster.yaml` values (picard/GATK 96 GB, fgbio 64 GB, bwa 16 CPUs / 36 GB); `conf/test.config` caps everything at 4 CPUs / 12 GB.
- **Session 2 starting point:** the `consensus` subworkflow can call `BAM_ALIGN` with the consensus unmapped BAM and `ext.prefix = "${meta.id}_consensus"`; install `fgbio/groupreadsbyumi`, `fgbio/callduplexconsensusreads`, `fgbio/filterconsensusreads`, `picard/markduplicates`, `samtools/stats`, `samtools/flagstat`, `qualimap/bamqc` and `multiqc` with the same `nf-core modules install` recipe (needs the temporary `.nf-core.yml`).
