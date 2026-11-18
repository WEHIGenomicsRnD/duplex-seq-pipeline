#!/usr/bin/env bash
#
# Build a one-genome FastQ Screen database (bowtie2 index + config) from a FASTA file.
#
# Usage: build_fastq_screen_db.sh <fasta[.gz]> <outdir> [name]
#
# Creates <outdir>/FastQ_Screen_Genomes/ with the bowtie2 index and a fastq_screen.conf whose
# DATABASE path is relative to that folder name, so the folder can be passed as
# `--fastq_screen_db <outdir>/FastQ_Screen_Genomes` and staged into a task directory.
#
# Only needs `bowtie2-build` on PATH (for example from conda: `bioconda::bowtie2`). Set
# BOWTIE2_BUILD to use a different command, such as a container wrapper.

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $(basename "$0") <fasta[.gz]> <outdir> [name]" >&2
    exit 1
fi

fasta=$1
outdir=$2
name=${3:-ecoli_atcc10798}
bowtie2_build=${BOWTIE2_BUILD:-bowtie2-build}

if [[ ! -f "$fasta" ]]; then
    echo "FASTA file not found: $fasta" >&2
    exit 1
fi

if ! command -v "${bowtie2_build%% *}" > /dev/null 2>&1; then
    echo "'${bowtie2_build%% *}' not found on PATH; install bowtie2 (for example 'conda install -c bioconda bowtie2') or set BOWTIE2_BUILD" >&2
    exit 1
fi

db_dir="${outdir}/FastQ_Screen_Genomes"
mkdir -p "${db_dir}/${name}"

# bowtie2-build reads plain or gzipped FASTA; the index prefix is <name>/<name>
${bowtie2_build} --threads 2 "$fasta" "${db_dir}/${name}/${name}"

printf 'THREADS\t2\nDATABASE\t%s\tFastQ_Screen_Genomes/%s/%s\tBOWTIE2\n' \
    "$name" "$name" "$name" > "${db_dir}/fastq_screen.conf"

echo "FastQ Screen database written to ${db_dir}"
