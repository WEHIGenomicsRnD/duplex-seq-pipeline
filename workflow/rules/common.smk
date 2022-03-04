import os
from glob import iglob

fqs = iglob("fastq/*_R1.fastq.gz")

# extract basename of full file path
base = [os.path.basename(i) for i in fqs]

# extract sample name
samples = []
for f in base:
    samples.append(f.split("_R1")[0])

correct_umis = config["umis"] != ""
umi_postfix = "uncorrected_" if correct_umis else ""


def get_create_dict_output():
    basename = os.path.splitext(config["ref"])[0]
    output = f"{basename}.dict"
    return output


def get_varcall_output():
    output = expand("results/variants/{sample}.vcf", sample=samples)
    return output


def get_calc_rinfo_output():
    output = expand("results/QC/read_info/{sample}.txt", sample=samples)
    return output
