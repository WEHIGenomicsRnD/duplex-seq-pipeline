#!/usr/bin/env python3
"""
Extract per-read duplex sequencing info from a BAM file.

Outputs tab-separated: query_name, chrom, pos, mpos, strand, umi
to stdout.
"""
import pysam
import os
import sys
from argparse import ArgumentParser

PROGRAM_NAME = 'get_read_info'

def parse_args():
    parser = ArgumentParser(description='Extract read info for duplex metrics.')
    parser.add_argument('bamfile', type=str, help='Input BAM file.')
    return parser.parse_args()

def get_read_info(bam):
    for read in bam:
        if read.reference_start < 0:
            continue
        if read.mate_is_unmapped:
            continue

        strand = ""
        if read.is_read1 and read.reference_start < read.next_reference_start:
            if not read.is_reverse and read.mate_is_reverse:
                strand = "+"
        elif read.is_read2 and read.reference_start > read.next_reference_start:
            if read.is_reverse and not read.mate_is_reverse:
                strand = "+"
        elif read.is_read1 and read.reference_start > read.next_reference_start:
            if read.is_reverse and not read.mate_is_reverse:
                strand = "-"
        elif read.is_read2 and read.reference_start < read.next_reference_start:
            if not read.is_reverse and read.mate_is_reverse:
                strand = "-"

        if strand not in ["-", "+"]:
            continue

        pos = min(read.reference_start, read.next_reference_start)
        mpos = max(read.reference_start, read.next_reference_start)
        umi = read.get_tag("ZB") + "-" + read.get_tag("ZA") if strand == "-" else read.get_tag("RX")

        print(f'{read.query_name}\t{read.reference_name}\t{pos}\t{mpos}\t{strand}\t{umi}',
              file=sys.stdout)

def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bamfile)
    sample = os.path.basename(args.bamfile)
    print("Processing", sample, "...", file=sys.stderr)
    get_read_info(bam)

if __name__ == "__main__":
    main()
