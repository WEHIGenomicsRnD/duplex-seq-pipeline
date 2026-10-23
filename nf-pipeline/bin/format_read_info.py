#!/usr/bin/env python3
"""
Reformat read info output to resemble the Sanger NanoSeq pipeline format.
Output: chrom, pos, mpos, umi, x, y (+ and - strand counts).

CLI replacement for the Snakemake-integrated script.
"""
import numpy as np
import pandas as pd
import sys
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser(description='Format read info into NanoSeq-like output.')
    parser.add_argument('input', type=str, help='Input read info summary file.')
    parser.add_argument('output', type=str, help='Output gzipped TSV file.')
    return parser.parse_args()

def main():
    args = parse_args()

    rinfo_summary = pd.read_csv(
        args.input,
        sep='\t',
        header=None,
        names=['chrom', 'start', 'end', 'strand', 'umi']
    )
    rinfo_summary['chrom'] = rinfo_summary.chrom.str.strip().values

    # re-split the output from uniq --count (leading/trailing spaces)
    tmp = rinfo_summary.chrom.str.split(pat=r'\s+', expand=True)
    rinfo_summary['chrom'] = tmp.loc[:, 1].values
    rinfo_summary['count'] = tmp.loc[:, 0].values
    del tmp

    # get both strands on the same row
    rinfo_summary = pd.merge(
        rinfo_summary[rinfo_summary.strand == '+'],
        rinfo_summary[rinfo_summary.strand == '-'],
        how='outer',
        on=['chrom', 'start', 'end', 'umi']
    ).drop_duplicates().fillna(0)

    rename_cols = {'count_x': 'x', 'count_y': 'y', 'start': 'pos', 'end': 'mpos'}
    rinfo_summary = rinfo_summary[['chrom', 'start', 'end', 'umi', 'count_x', 'count_y']]
    rinfo_summary = rinfo_summary.rename(mapper=rename_cols, axis=1)
    rinfo_summary.to_csv(args.output, sep='\t', index=False, compression='gzip')

if __name__ == "__main__":
    main()
