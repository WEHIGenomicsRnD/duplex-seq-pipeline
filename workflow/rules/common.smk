import os
from glob import iglob

fqs = iglob('fastq/*_R1.fastq.gz')

# extract basename of full file path
base = [os.path.basename(i) for i in fqs]

# extract sample name
samples = []
for f in base:
    samples.append(f.split('_R1')[0])

def get_filtered_consensus_output():
    output = expand(
        'results/consensus/{sample}_consensus_mapped_merged_filtered.bam',
        sample=samples
    )
    return output
