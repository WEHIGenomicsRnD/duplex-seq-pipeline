"""
Reformats read info output to be similar to the Sanger
NanoSeq pipeline (https://github.com/cancerit/NanoSeq)
The format is: chrom, pos, mpos, epos, umi, x, y
where x and y are read counts on the + and - strands
respectively
"""
from argparse import ArgumentParser

import pandas as pd


def parse_args():
    '''Parse arguments'''
    description = '''
        Reformats read info to create duplex family counts. Takes as input:
        id chrom pos mpos epos strand umi
        to:
        chrom pos mpos epos umi x y
        where x and y are read counts on the + and - strands respectively
        '''
    parser = ArgumentParser(description=description)
    parser.add_argument('input_file',
                        metavar='INPUT_FILE',
                        type=str,
                        help='Input rinfo_summary file.')
    parser.add_argument('output_file',
                        metavar='OUTPUT_FILE',
                        type=str,
                        help='Output final rinfo file.')
    return parser.parse_args()


def main(input_file, output_file):
    rinfo_summary = pd.read_csv(input_file,
                                sep='\t',
                                header=None,
                                names=['chrom', 'pos', 'mpos', 'epos', 'strand', 'umi'])
    rinfo_summary['chrom'] = rinfo_summary.chrom.str.strip().values

    # we have to resplit the output from uniq --count due to the leading and trailing
    # spaces surrounding the count number. We could do this above, with read.csv,
    # but it's slower as regex split forces the python engine to load the file.
    tmp = rinfo_summary.chrom.str.split(pat="\\s+", expand=True)
    rinfo_summary['chrom'] = tmp.loc[:, 1].values
    rinfo_summary['count'] = tmp.loc[:, 0].values
    del tmp  # del tmp from memory as it can be large

    # get both strands on the same row
    rinfo_summary = pd.merge(
        rinfo_summary[rinfo_summary.strand == '+'],
        rinfo_summary[rinfo_summary.strand == '-'],
        how='outer',
        on=['chrom', 'pos', 'mpos', 'epos', 'umi']
    ).drop_duplicates().fillna(0)

    # some final cleaning and output to file
    rename_cols = {'count_x': 'x', 'count_y': 'y'}
    rinfo_summary = rinfo_summary[['chrom', 'pos', 'mpos', 'epos', 'umi', 'count_x', 'count_y']]
    rinfo_summary = rinfo_summary.rename(mapper=rename_cols, axis=1)
    rinfo_summary.to_csv(output_file, sep='\t', index=False, compression='gzip')


if __name__ == "__main__":
    args = parse_args()
    input_file = args.input_file
    output_file = args.output_file

    main(input_file, output_file)
