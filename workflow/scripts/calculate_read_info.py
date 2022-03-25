"""
Script for generating stats on duplex reads.
Counts the number of families, their min, mean,
media, and max, as well as the proportion of
families with reads > 1, the number of paired
families (+ and - strand) and number of paired
families where each family has reads > 1.
"""
import pysam
import numpy as np
import pandas as pd
import seaborn as sns
import os
import sys

sys.stderr = open(snakemake.log[0], "w")

bamfile = snakemake.input[0]
outfile = snakemake.output[0]

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

        r = (read.query_name,
             read.reference_name,
             pos,
             mpos,
             strand,
             umi)
        rinfo.append(r)

        if len(rinfo) % 1000000 == 0:
            print("Processed %d reads" % len(rinfo), file=sys.stdout)

    rinfo = pd.DataFrame(rinfo, columns=["rname", "chrom", "pos", "mpos", "strand", "umi"]).drop_duplicates()
    return rinfo

def process_reads(bamfile):

    bam = pysam.AlignmentFile(bamfile)

    sample = os.path.basename(bamfile)
    print("Processing", sample, "...", file=sys.stdout)

    rinfo = get_read_info(bam)
    rinfo_summary = rinfo.groupby(["chrom", "pos", "mpos", "strand", "umi"]).size().reset_index()

    return rinfo_summary

def main():
    rinfo_summary = process_reads(bamfile)

    # restructure format to resemble the Sanger NanoSeq pipeline output
    rinfo_summary = pd.merge(
        rinfo_summary[rinfo_summary.strand == '+'],
        rinfo_summary[rinfo_summary.strand == '-'],
        how='outer',
        on=['chrom', 'pos', 'mpos', 'umi']
    ).drop_duplicates().fillna(0)

    rinfo_summary = rinfo_summary[['chrom', 'pos', 'mpos', 'umi', '0_x', '0_y']]
    rinfo_summary = rinfo_summary.rename(mapper={'0_x': 'x', '0_y': 'y'}, axis=1)

    rinfo_summary.to_csv(outfile, sep="\t", index=False, compression="gzip")

if __name__ == "__main__":
    main()
