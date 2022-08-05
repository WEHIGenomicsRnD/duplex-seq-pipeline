"""
Reformats read info output to be similar to the Sanger
NanoSeq pipeline (https://github.com/cancerit/NanoSeq)
The format is: chrom, start, end, UMI, x, y
where x and y are read counts on the + and - strands
respectively
"""
import numpy as np
import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

rinfo_summary_file = snakemake.input[0]
outfile = snakemake.output[0]

def main():
    rinfo_summary = pd.read_csv(rinfo_summary_file,
                                sep='\t',
                                header=None,
                                names=['chrom', 'start', 'end', 'strand', 'umi'])
    rinfo_summary['chrom'] = rinfo_summary.chrom.str.strip().values

    # we have to resplit the output from uniq --count due to the leading and trailing
    # spaces surrounding the count number. We could do this above, with read.csv,
    # but it's slower as regex split forces the python engine to load the file.
    tmp = rinfo_summary.chrom.str.split(pat='\s+', expand=True)
    rinfo_summary['chrom'] = tmp.loc[:,1].values
    rinfo_summary['count'] = tmp.loc[:,0].values
    del tmp # del tmp from memory as it can be large

    # get both strands on the same row
    rinfo_summary = pd.merge(
        rinfo_summary[rinfo_summary.strand == '+'],
        rinfo_summary[rinfo_summary.strand == '-'],
        how='outer',
        on=['chrom', 'start', 'end', 'umi']
    ).drop_duplicates().fillna(0)

    # some final cleaning and output to file
    rename_cols = {'count_x': 'x', 'count_y': 'y', 'start': 'pos', 'end': 'mpos'}
    rinfo_summary = rinfo_summary[['chrom', 'start', 'end', 'umi', 'count_x', 'count_y']]
    rinfo_summary = rinfo_summary.rename(mapper=rename_cols, axis=1)
    rinfo_summary.to_csv(outfile, sep='\t', index=False, compression='gzip')

if __name__ == "__main__":
    main()
