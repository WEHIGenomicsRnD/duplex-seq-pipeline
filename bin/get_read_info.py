"""
Script for generating stats on duplex reads.
Outputs each read's query name, reference it's
aligned to, position, the mate read's end position,
as well as the strand-adjusted UMI (the UMI order
is always in the forward orientation, even if on
the reverse strand).
"""
import pysam
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
    mate_cache = {}
    for read in bam:

        # Skip unmapped and unpaired reads
        if read.is_unmapped or read.mate_is_unmapped or not read.is_paired:
            continue

        # Skip secondary and supplementary alignments
        if read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue

        # Filter out reads that are not on the same chromosome
        if read.reference_id != read.next_reference_id:
            continue

        qn = read.query_name
        prev = mate_cache.pop(qn, None)

        if prev is None:
            mate_cache[qn] = read
            continue

        r1, r2 = (read, prev) if read.is_read1 else (prev, read)

        left_is_r1 = r1.reference_start <= r2.reference_start
        if left_is_r1:
            # FR expected: r1 forward, r2 reverse
            strand = "+" if (not r1.is_reverse and r2.is_reverse) else ""
        else:
            # RF expected: r1 reverse, r2 forward
            strand = "-" if (r1.is_reverse and not r2.is_reverse) else ""

        if not strand:
            continue

        if strand == "-":
            if not (r1.has_tag("ZA") and r1.has_tag("ZB")):
                continue
            umi = f'{r1.get_tag("ZB")}-{r1.get_tag("ZA")}'
        else:
            if not r1.has_tag("RX"):
                continue
            umi = r1.get_tag("RX")

        pos = min(r1.reference_start, r2.reference_start)
        mpos = max(r1.reference_start, r2.reference_start)
        epos = max(r1.reference_end, r2.reference_end)

        print(f'{r1.query_name}\t{r1.reference_name}\t{pos}\t{mpos}\t{epos}\t{strand}\t{umi}',
              file=sys.stdout)


def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.bamfile)

    sample = os.path.basename(args.bamfile)
    print("Processing", sample, "...", file=sys.stderr)

    get_read_info(bam)


if __name__ == "__main__":
    main()
