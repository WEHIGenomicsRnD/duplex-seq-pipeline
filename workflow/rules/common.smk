import os
from glob import iglob

fqs = iglob('fastq/*_R1.fastq.gz')

# extract basename of full file path
base = [os.path.basename(i) for i in fqs]

# extract sample name
samples = []
for f in base:
    samples.append(f.split('_R1')[0])

def get_varcall_output():
    output = expand(
        'results/variants/{sample}.vcf',
        sample=samples
    )
    return output

def get_qualimap_output():
    output = expand(
        'results/QC/qualimap/{sample}/qualimapReport.html',
        sample=samples
    )
    return output

def get_assembly_output():
    output = expand(
        'results/QC/checkm/assembly/{sample}.contig',
        sample=samples
    )
    return output
