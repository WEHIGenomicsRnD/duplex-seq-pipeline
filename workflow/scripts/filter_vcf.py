"""
Add DP field in the vcf file
"""
#import numpy as np
#import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

vcf_file = snakemake.input[0]
outfile = snakemake.output[0]


def main():
    ofile=open(outfile,'w')
    flag=1
    with open (vcf_file) as vcf:
       for l in vcf.readlines():
          if l.startswith("#"):
              if flag and l.startswith("##FORMAT"):
                  ofile.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">' + '\n')
                  flag=0
              ofile.write(l)
          else:
              if l.strip().split("\t")[4] == "N": continue     ## Skipping rows with alt allele N
              dp_res=count_depth(l)
              ofile.write("\t".join(dp_res) + '\n')



def count_depth(l):
    line=l.strip().split("\t")

    if "N" in line[4].strip().split(",") and len(line[4].split(",")) > 1:
        line=process_allele(l)

    nc=l.strip().split("\t")[-1].split(":")[-1].split(",")

    dp=0
    for row in nc:
        if ("d" not in row) and (row != "") and (row != "."):
            dp=int(row.split("=")[-1]) + dp
    line[7]=line[7]+";DP="+str(dp)
    line[9]=line[9]+":"+str(dp)
    line[8]=line[8]+":"+"DP"

    return line




def process_allele(line):
    row=line.strip().split("\t")
    allele_base=row[4].strip().split(',')
    afreq=row[7].strip().split(";")[-1].split("=")[-1].split(",")
    ac=row[7].strip().split(";")[0].split("=")[-1].split(",")
#    print(f"AC - {ac}")
    for index, base in enumerate(allele_base):
        if base =="N":
           rm_base=allele_base.pop(index)
           rm_freq=afreq.pop(index)
           rm_ac=ac.pop(index)

           new_afreq="AF="+",".join(afreq)
           new_ac="AC="+",".join(ac)
           row[7]=";".join([new_ac,new_afreq])
           row[4]=",".join(allele_base)
    return row


if __name__ == "__main__":
    main()
