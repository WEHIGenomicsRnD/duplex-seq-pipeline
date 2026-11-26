#!/usr/bin/env python3
"""
Turn the bcftools mpileup records of the candidate variant positions into the VCF that the
naive variant caller (nvc) followed by filter_vcf.py writes.

bcftools mpileup writes one record for the nucleotides of a position and, where reads have
gaps, a second record for the indels anchored on it. The records of a position are merged
into one record per position, with the alleles named and ordered as nvc does:

    * the nucleotide counts are the AD values of the nucleotide record, and the N count
        (bases called N, which mpileup leaves out of AD) is its DP minus the sum of AD,
    * an insertion is counted under the anchor base followed by the inserted bases, a
        deletion under its length,
    * alleles (the reference allele included) supported by fewer than --min_support_depth
        reads are dropped; the remaining alleles other than the reference are the ALT alleles,
        most supported first,
    * REF is extended over the longest deletion, and the other alleles with it.

The INFO and FORMAT fields are those of nvc (AC, AF; GT, AC, AF, NC), and as filter_vcf.py
does, N alleles are removed and the read depth (DP) is added. AF is the allele count over the
summed counts of the alleles that passed --min_support_depth, as in nvc.
"""

import gzip
from argparse import ArgumentParser
from itertools import groupby

# 0-based indices of the VCF columns used below
CHROM_COLUMN = 0
POS_COLUMN = 1
REF_COLUMN = 3
ALT_COLUMN = 4
INFO_COLUMN = 7
FORMAT_COLUMN = 8
SAMPLE_COLUMN = 9

# mpileup's allele for "any other allele", which never has reads
UNOBSERVED_ALLELE = "<*>"
UNKNOWN_NUCLEOTIDE = "N"

HEADER_LINES = [
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">\n',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
    '##FORMAT=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n',
    '##FORMAT=<ID=AF,Number=.,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">\n',
    '##FORMAT=<ID=NC,Number=.,Type=String,Description="Nucleotide and indel counts">\n',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n',
]


def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(
        description="Merge the bcftools mpileup records of each position into an nvc-style VCF."
    )
    parser.add_argument(
        "input_vcf",
        metavar="INPUT_VCF",
        help="VCF (optionally gzipped) written by bcftools mpileup with FORMAT/AD and FORMAT/DP.",
    )
    parser.add_argument("output_vcf", metavar="OUTPUT_VCF", help="Filtered VCF.")
    parser.add_argument(
        "--min_support_depth",
        type=int,
        default=2,
        help="Minimum number of reads supporting an allele (default: 2).",
    )
    parser.add_argument(
        "--ploidy", type=int, default=1, help="Ploidy of the genotype (default: 1)."
    )
    return parser.parse_args()


def format_float(value):
    """Format a float as Python 2's str() did (12 significant digits), as nvc writes it."""
    text = f"{value:.12g}"
    if not any(c in text for c in ".einf"):
        text += ".0"
    return text


def sort_key(allele_count):
    """Sort key of nvc for (allele, count): by count, then by allele name.

    Under Python 2, which nvc runs on, a deletion (an int) sorts before any string, so with
    reverse=True the deletions come after the nucleotides and insertions of the same count.
    """
    allele, count = allele_count
    return count, isinstance(allele, str), allele


def sample_values(fields):
    """Return the FORMAT values of the first sample as a dict."""
    keys = fields[FORMAT_COLUMN].split(":")
    values = fields[SAMPLE_COLUMN].split(":")
    return dict(zip(keys, values))


def trim_common_suffix(ref, alt):
    """Remove the bases shared at the end of REF and ALT, keeping at least one base of each."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    return ref, alt


def count_alleles(records):
    """Count the reads of each allele at a position from its mpileup records.

    Returns the reference base, the reference sequence covering the longest deletion seen and
    a dict of read counts keyed as nvc does: a nucleotide by its base, an insertion by the
    anchor base followed by the inserted bases, a deletion by its length (int).
    """
    ref_base = None
    ref_sequence = ""
    counts = {}
    for fields in records:
        # nvc reads the reference in upper case
        ref = fields[REF_COLUMN].upper()
        alts = [a if a == UNOBSERVED_ALLELE else a.upper() for a in fields[ALT_COLUMN].split(",")]
        values = sample_values(fields)
        allele_depths = [int(x) for x in values["AD"].split(",")]
        ref_base = ref[0]
        if len(ref) > len(ref_sequence):
            ref_sequence = ref

        is_indel = "INDEL" in fields[INFO_COLUMN].split(";")
        if not is_indel:
            counts[ref] = counts.get(ref, 0) + allele_depths[0]
            for alt, depth in zip(alts, allele_depths[1:]):
                if alt != UNOBSERVED_ALLELE and depth:
                    counts[alt] = counts.get(alt, 0) + depth
            # bases called N are in DP but not in AD
            unknown = int(values["DP"]) - sum(allele_depths)
            if unknown > 0:
                counts[UNKNOWN_NUCLEOTIDE] = counts.get(UNKNOWN_NUCLEOTIDE, 0) + unknown
            continue

        for alt, depth in zip(alts, allele_depths[1:]):
            if alt == UNOBSERVED_ALLELE or not depth:
                continue
            trimmed_ref, trimmed_alt = trim_common_suffix(ref, alt)
            if len(trimmed_alt) > len(trimmed_ref):
                name = trimmed_alt[: len(trimmed_alt) - len(trimmed_ref) + 1]
            else:
                name = len(trimmed_ref) - len(trimmed_alt)
            counts[name] = counts.get(name, 0) + depth
    return ref_base, ref_sequence, counts


def allele_sequence(name, ref_base, indeled_ref):
    """The sequence of an allele on the extended reference, as nvc's _rework_indels builds it."""
    if isinstance(name, int):
        return indeled_ref[:1] + indeled_ref[1 + name :]
    if len(name) > 1:
        return indeled_ref[:1] + name[1:] + indeled_ref[1:]
    return name + indeled_ref[1:]


def call_genotype(passing, ploidy):
    """The genotype calls of nvc: the alleles with the most reads, one per ploidy."""
    if not passing:
        return [None] * ploidy
    nucs_sum = sum(count for _, count in passing)
    count_per_allele = float(nucs_sum) / ploidy
    nucs = [(count / count_per_allele, sequence) for sequence, count in passing]
    calls = []
    for _ in range(ploidy):
        nucs = sorted(nucs, reverse=True)
        for i, (weight, sequence) in enumerate(nucs):
            j = i + 1
            if sequence not in calls or j >= len(nucs) or weight > nucs[j][0]:
                break
        nucs[i] = (nucs[i][0] - 1, nucs[i][1])
        calls.append(sequence)
    return calls


def nc_entry(name, count):
    """An entry of nvc's NC field: deletions are written as d<length>."""
    label = f"d{name}" if isinstance(name, int) else name
    return f"{label}={count}"


def build_record(chrom, pos, records, min_support_depth, ploidy):
    """Return the nvc-style VCF line of a position, or None if it has no ALT allele."""
    ref_base, ref_sequence, counts = count_alleles(records)

    passing = {name: count for name, count in counts.items() if count >= min_support_depth}
    max_deletion = max((n for n in passing if isinstance(n, int)), default=0)
    indeled_ref = ref_base + ref_sequence[1 : 1 + max_deletion]

    # ALT alleles: the passing alleles other than the reference base, most supported first;
    # N alleles are left out (as filter_vcf.py does)
    ordered = sorted(passing.items(), key=sort_key, reverse=True)
    alts = []
    alt_counts = []
    for name, count in ordered:
        if name == ref_base or name == UNKNOWN_NUCLEOTIDE:
            continue
        alts.append(allele_sequence(name, ref_base, indeled_ref))
        alt_counts.append(count)
    if not alts:
        return None

    # AF is over all passing alleles: the reference and N included
    nucs_sum = sum(passing.values())
    allele_freqs = [format_float(float(count) / nucs_sum) for count in alt_counts]

    passing_sequences = [
        (allele_sequence(name, ref_base, indeled_ref), count) for name, count in ordered
    ]
    alleles = [indeled_ref] + alts
    genotype = "/".join(
        str(alleles.index(call)) if call in alleles else "."
        for call in call_genotype(passing_sequences, ploidy)
    )

    # read depth as filter_vcf.py sums it from NC: all counts but the deletions
    depth = sum(count for name, count in counts.items() if not isinstance(name, int))
    nc = "".join(
        nc_entry(name, count) + ","
        for name, count in sorted(counts.items(), key=lambda x: str(x[0]))
        if count
    )

    ac = ",".join(str(c) for c in alt_counts)
    af = ",".join(allele_freqs)
    info = f"AC={ac};AF={af};DP={depth}"
    sample = f"{genotype}:{ac}:{af}:{nc}:{depth}"
    return "\t".join(
        [chrom, str(pos), ".", indeled_ref, ",".join(alts), ".", ".", info, "GT:AC:AF:NC:DP", sample]
    )


def open_vcf(path):
    """Open a plain or gzipped VCF for reading."""
    with open(path, "rb") as handle:
        is_gzipped = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gzipped else open(path)


def filter_pileup(input_vcf, output_vcf, min_support_depth, ploidy):
    """Write the nvc-style VCF of the mpileup records in input_vcf."""
    with open_vcf(input_vcf) as vcf_in, open(output_vcf, "w") as vcf_out:
        # header: the contig lines of mpileup, then the nvc fields and DP
        vcf_out.write("##fileformat=VCFv4.3\n##source=filter_pileup.py (bcftools mpileup)\n")
        for line in vcf_in:
            if line.startswith(("##contig", "##reference")):
                vcf_out.write(line)
            elif line.startswith("#CHROM"):
                vcf_out.writelines(HEADER_LINES)
                vcf_out.write(line)
                break

        # mpileup writes the records of a position next to each other
        records = (line.rstrip("\n").split("\t") for line in vcf_in if line.strip())
        for (chrom, pos), position_records in groupby(
            records, key=lambda f: (f[CHROM_COLUMN], int(f[POS_COLUMN]))
        ):
            line = build_record(chrom, pos, list(position_records), min_support_depth, ploidy)
            if line is not None:
                vcf_out.write(line + "\n")


def main():
    args = parse_args()
    filter_pileup(args.input_vcf, args.output_vcf, args.min_support_depth, args.ploidy)


if __name__ == "__main__":
    main()
