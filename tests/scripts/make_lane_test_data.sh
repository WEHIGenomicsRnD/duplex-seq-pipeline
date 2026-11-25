#!/usr/bin/env bash
#
# Cut the lane test data in tests/data/lanes from the test FASTQ files in tests/data.
#
# Usage: make_lane_test_data.sh [fastq_dir]
#
# NanoMB1Rep1 gets three lanes (L001 to L003) of 3,000 read pairs each, and NanoMB1Rep2 one lane (L001) of
# 3,000 read pairs, so that the lanes of a sample are consecutive pieces of the same file and the result of
# merging them is easy to predict. The files are written to <fastq_dir>/lanes.

set -euo pipefail

fastq_dir=${1:-"$(dirname "$0")/../data"}
source_name=NanoMB1Rep1_HJK2GDSX3_sampled
pairs_per_lane=3000
lines_per_lane=$((pairs_per_lane * 4))

mkdir -p "${fastq_dir}/lanes"

# cut_lane <sample> <lane number> <first pair, counted from 0>
cut_lane() {
    local sample=$1 lane=$2 first_pair=$3
    local first_line=$((first_pair * 4 + 1))
    local last_line=$((first_line + lines_per_lane - 1))
    for read in R1 R2; do
        # gzip -n keeps no timestamp, so the files are identical every time
        zcat "${fastq_dir}/${source_name}_${read}.fastq.gz" \
            | sed -n "${first_line},${last_line}p" \
            | gzip -n -9 > "${fastq_dir}/lanes/${sample}_L$(printf '%03d' "${lane}")_${read}.fastq.gz"
    done
}

cut_lane NanoMB1Rep1 1 0
cut_lane NanoMB1Rep1 2 "${pairs_per_lane}"
cut_lane NanoMB1Rep1 3 $((2 * pairs_per_lane))
cut_lane NanoMB1Rep2 1 $((3 * pairs_per_lane))
