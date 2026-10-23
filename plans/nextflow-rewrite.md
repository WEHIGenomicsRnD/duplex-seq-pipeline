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

| Commit | Change |
|---|---|
| 2960c23 | `run_qualimap` param (default false); qualimap becomes optional |
| c1550c2 | `get_read_info.py` rewritten: mate pairing, stricter filters, `epos` column, `sort -k 2,7` / `cut -f 2-7` |
| e0a2a3f | read info merged into one step, with an empty-output check |
| 277c506 | `target_bed` param, passed as `--target_bed` to calc-duplex-metrics |
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
- **Versions:** upgrade to current nf-core module versions (fgbio 2.5.x, picard 3.x, etc.). Validate against Snakemake on the test data and record the differences in CHANGELOG/docs.
  - Check after the upgrade: fgbio 2.x GroupReadsByUmi still accepts a coordinate-sorted MBA BAM (it sorts internally if needed).
  - Check after the upgrade: `ExtractUmisFromBam` still exists and keeps `--molecular-index-tags`.

## Containers needed

| Tool | Container source |
|---|---|
| fastqc, fastq-screen (+bowtie2), picard, fgbio, bwa+samtools, samtools, qualimap | nf-core modules (Seqera/biocontainers); pinned biocontainer tags if we keep the Snakemake versions |
| varscan + samtools, nvc + samtools | multi-tool: Seqera Containers built from the bioconda spec (nvc 0.0.4 is on bioconda) |
| python/pysam/pandas (`get_read_info.py`, `format_read_info.py`, `filter_vcf.py`) | Seqera Containers built from the conda spec |
| calc-duplex-metrics | **exists**: `ghcr.io/wehigenomicsrnd/calculate-duplex-metrics:<release tag>` (build-docker.yml in that repo). Built on release. Latest is `0.2.0`; check that it includes BAM mode and 39fc6c9, otherwise cut `0.2.1`. `calc-duplex-metrics` is on PATH at `/usr/local/bin`. No runtime install. |
| MultiQC + duplex-multiqc plugin | **needs building**: add a Dockerfile (`FROM` the Seqera multiqc 1.3x image, `pip install duplex-multiqc@<tag>`) and a GHCR build workflow in the duplex-multiqc repo. Conda profile: `multiqc` + pip install from git. |

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
  - MergeBamAlignment `SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR CREATE_INDEX=true`
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

| Area | Modules |
|---|---|
| QC | `fastqc`, `fastqscreen/fastqscreen`, `qualimap/bamqc`, `multiqc` |
| Picard | `picard/createsequencedictionary`, `picard/fastqtosam`, `picard/mergebamalignment`, `picard/markduplicates` |
| Alignment and samtools | `bwa/index`, `bwa/mem` (sort=false, `-p`), `samtools/faidx`, `samtools/stats`, `samtools/flagstat`, `samtools/index` |
| fgbio | `fgbio/groupreadsbyumi`, `fgbio/callduplexconsensusreads`, `fgbio/filterconsensusreads` |

- **Wiring:** check each signature after install and pass params through `val` inputs or `ext.args`, so that every Snakemake argument reaches its tool (this fixes the skeleton's mismatched calls).
- **Local modules** (meta map, `versions.yml` or topic `versions`, `stub:` block, `ext.args`/`ext.prefix`, no `params.*` inside the module):
  - **`FGBIO_EXTRACTUMISFROMBAM`, `FGBIO_CORRECTUMIS`, `FGBIO_CLIPBAM`:** Seqera fgbio container (or an nf-core module if one exists).
  - **`PICARD_SAMTOFASTQ`:** local only if no nf-core module exists; picard biocontainer, `INTERLEAVE=true`.
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
- **Why not fastquorum's data:** it uses `10M1S+T` and chr17, which doesn't match the `3M2S146T` structure, and its samples aren't usable for duplex metrics here.
- **Source:** the user's ~30 GB E. coli duplex sample, which isn't uploaded anywhere yet.
- **Approach: take a region, not a random fraction.**
  - Random `samtools view -s` thins every UMI family binomially. With `min_reads '4 2 2'` (at least 2 reads per strand), a 0.1% sample gives almost no duplexes, so the duplex metrics and VCF would be empty.
  - Taking every read pair from a small genomic window keeps the families intact at their full depth.
- **Steps:**
  1. From the Snakemake `results/mapped_bams/{s}_mapped_merged.bam` (coordinate-sorted, indexed), pick a window big enough for roughly **30–50k read pairs**:
     - window size ≈ target_pairs ÷ (total_pairs ÷ 4.6 Mb);
     - for example, 150M pairs gives about 33k pairs/kb, so a window of about **1–2 kb**.

     Then run `samtools view -F 0x900 bam chr:start-end | cut -f1 | sort -u > names.txt`.
  2. The mapped BAM has the UMI bases removed, so take those names from the **original FASTQs** with `seqkit grep -f names.txt` (or `samtools view -N names.txt` on the unmapped BAM, followed by `samtools fastq`). That keeps the `3M2S` UMI bases. Output: about 2 × 3–5 MB gzipped FASTQ.
  3. Reference: the full BL21 genome (about 4.6 MB, about 1.4 MB gzipped). bwa index is fast, and off-target mates map properly.
  4. `target.bed`: a sub-window inside the extracted region (exercises `--target_bed`). `umis.txt`: the existing whitelist (`metadata/umis.txt`), if there is one, to exercise CorrectUmis.
  5. Optionally a second small sample (a different window, or a 2-sample samplesheet) to test multi-sample MultiQC merging.
- **Hosting:** the user's existing repo. `conf/test.config` points at raw URLs there. Target: `-profile test,docker` finishes in under 5 minutes on a GitHub runner.

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
3. **Real run:** `nextflow run . -profile test,docker` on the downsampled E. coli data. Every step should run and a single `multiqc/multiqc_report.html` should be produced, with R1/R2/raw/consensus grouped per sample in General Statistics and the duplex metrics section present.
4. **Parity with Snakemake:** run the Snakemake pipeline on the same test data and compare:
   - consensus BAM read counts (`samtools flagstat`)
   - `*_metrics.csv` (should be identical)
   - filtered VCF records
   - read-info output with `--output_read_info`
5. **Optional branches:** `--run_qualimap`, `--output_read_info`, `--umis`, `--target_bed`, `--variant_caller varscan`.
6. `nf-test test --profile docker` passes locally and in CI.
