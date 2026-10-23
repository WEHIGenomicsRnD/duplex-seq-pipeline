#!/usr/bin/env python3
"""
Add DP (read depth) field to a VCF file and filter out rows
with alt allele 'N'.

CLI replacement for the Snakemake-integrated script.
"""
import sys
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser(description='Add DP field and filter VCF.')
    parser.add_argument('input', type=str, help='Input VCF file.')
    parser.add_argument('output', type=str, help='Output filtered VCF file.')
    return parser.parse_args()

def count_depth(line):
    cols = line.strip().split("\t")

    if "N" in cols[4].strip().split(",") and len(cols[4].split(",")) > 1:
        cols = process_allele(line)

    nc = line.strip().split("\t")[-1].split(":")[-1].split(",")

    dp = 0
    for row in nc:
        if ("d" not in row) and (row != "") and (row != "."):
            dp = int(row.split("=")[-1]) + dp
    cols[7] = cols[7] + ";DP=" + str(dp)
    cols[9] = cols[9] + ":" + str(dp)
    cols[8] = cols[8] + ":" + "DP"

    return cols

def process_allele(line):
    row = line.strip().split("\t")
    allele_base = row[4].strip().split(',')
    afreq = row[7].strip().split(";")[-1].split("=")[-1].split(",")
    ac = row[7].strip().split(";")[0].split("=")[-1].split(",")

    for index, base in enumerate(allele_base):
        if base == "N":
            allele_base.pop(index)
            afreq.pop(index)
            ac.pop(index)

            new_afreq = "AF=" + ",".join(afreq)
            new_ac = "AC=" + ",".join(ac)
            row[7] = ";".join([new_ac, new_afreq])
            row[4] = ",".join(allele_base)
    return row

def main():
    args = parse_args()

    with open(args.input) as vcf, open(args.output, 'w') as ofile:
        flag = True
        for line in vcf:
            if line.startswith("#"):
                if flag and line.startswith("##FORMAT"):
                    ofile.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
                    flag = False
                ofile.write(line)
            else:
                if line.strip().split("\t")[4] == "N":
                    continue
                dp_res = count_depth(line)
                ofile.write("\t".join(dp_res) + '\n')

if __name__ == "__main__":
    main()
