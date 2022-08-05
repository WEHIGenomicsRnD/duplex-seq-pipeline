"""
Script for generating stats on duplex reads.
Outputs each read's query name, reference it's
aligned to, position, the mate read's end position,
as well as the strand-adjusted UMI (the UMI order
is always in the forward orientation, even if on
the reverse strand).
"""
import pysam
import numpy as np
import pandas as pd
import os
import sys
from argparse import ArgumentParser

PROGRAM_NAME = 'get_read_info'
EXIT_FILE_IO_ERROR = 1
EXIT_INPUT_ERROR = 2

def exit_with_error(message, exit_status):
    '''
    Print message and exit with exit_status
    '''
    print(f"{PROGRAM_NAME} ERROR: {message}, exiting", file=sys.stderr)
    sys.exit(exit_status)

def parse_args():
    '''Parse arguments'''
    description = '''
        Extract the following read info: name, chromosome, start,
        end, strand and UMI.
        '''
    parser = ArgumentParser(description=description)
    parser.add_argument('bamfile',
                        metavar='BAMFILE',
                        type=str,
                        help='Input BAM file.')
    return parser.parse_args()

def get_read_info(bam):
    rinfo = []
    for read in bam:

        if read.reference_start < 0:
            # read not mapped
            continue
        if read.mate_is_unmapped:
            continue

        strand = ""
        if read.is_read1 and read.reference_start < read.next_reference_start:
            # read 1 is left-most read
            if not read.is_reverse and read.mate_is_reverse:
                # check that orientation is proper
                strand = "+"
        elif read.is_read2 and read.reference_start > read.next_reference_start:
            # read 2 is right-most read
            if read.is_reverse and not read.mate_is_reverse:
                # same as above, we've just selected read 2
                strand = "+"
        elif read.is_read1 and read.reference_start > read.next_reference_start:
            # read 1 is right-most read
            if read.is_reverse and not read.mate_is_reverse:
                strand = "-"
        elif read.is_read2 and read.reference_start < read.next_reference_start:
            # read 2 is left-most read
            if not read.is_reverse and read.mate_is_reverse:
                # same as above, we've just selected read 2
                strand = "-"

        if strand not in ["-", "+"]:
            continue

        pos = min(read.reference_start, read.next_reference_start)
        mpos = max(read.reference_start, read.next_reference_start)
        umi = read.get_tag("ZB") + "-" + read.get_tag("ZA") if strand == "-" else read.get_tag("RX")

        print(f'{read.query_name}\t{read.reference_name}\t{pos}\t{mpos}\t{strand}\t{umi}', file=sys.stdout)

def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bamfile)

    sample = os.path.basename(args.bamfile)
    print("Processing", sample, "...", file=sys.stderr)

    get_read_info(bam)

if __name__ == "__main__":
    main()
