#!/usr/bin/env python3
"""
Post-process the VCF written by the naive variant caller (nvc).

For every variant record this script:

    * skips the record if its only ALT allele is "N" (there is no real variant),
    * removes an "N" allele from a multi-allelic record, together with its values
        in the per-allele AC and AF fields of the INFO column,
    * adds the read depth (DP) to the INFO column and to the sample column.

A FORMAT header line declaring DP is added to the header. The read depth is the sum
of the per-nucleotide counts that nvc writes in the last field of the sample column;
the deletion counts are not included. Only the first sample is annotated, as nvc
writes a single sample.
"""

from argparse import ArgumentParser

# 0-based indices of the VCF columns used below
ALT_COLUMN = 4
INFO_COLUMN = 7
FORMAT_COLUMN = 8
SAMPLE_COLUMN = 9

# INFO fields with one value per ALT allele
PER_ALLELE_INFO_KEYS = ("AC", "AF")

DP_HEADER_LINE = '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'


def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(
        description="Add the read depth (DP) to an nvc VCF and drop N alleles."
    )
    parser.add_argument("input_vcf", metavar="INPUT_VCF", help="VCF written by nvc.")
    parser.add_argument(
        "output_vcf", metavar="OUTPUT_VCF", help="Filtered VCF with DP added."
    )
    return parser.parse_args()


def sum_read_depth(sample_column):
    """Sum the per-nucleotide read counts in a sample column.

    nvc writes the counts in the last ":"-separated field of the sample column, as a
    comma-separated list of "<label>=<count>" entries. Entries containing a "d"
    (the deletion counts) are left out, as are empty and missing (".") entries.
    """
    count_entries = sample_column.split(":")[-1].split(",")

    depth = 0
    for entry in count_entries:
        if "d" in entry or entry in ("", "."):
            continue
        depth += int(entry.split("=")[-1])
    return depth


def remove_n_alleles(fields):
    """Remove "N" from a multi-allelic ALT column, in place.

    The matching values are also removed from the per-allele INFO fields (AC and AF);
    all other INFO fields are left as they are.
    """
    alt_alleles = fields[ALT_COLUMN].split(",")
    if len(alt_alleles) < 2 or "N" not in alt_alleles:
        return

    # positions of the alleles that are kept
    kept = [i for i, allele in enumerate(alt_alleles) if allele != "N"]
    fields[ALT_COLUMN] = ",".join(alt_alleles[i] for i in kept)

    info_fields = []
    for info_field in fields[INFO_COLUMN].split(";"):
        key, _, value = info_field.partition("=")
        if key in PER_ALLELE_INFO_KEYS:
            kept_values = [value.split(",")[i] for i in kept]
            info_field = f"{key}={','.join(kept_values)}"
        info_fields.append(info_field)
    fields[INFO_COLUMN] = ";".join(info_fields)


def add_read_depth(fields):
    """Add the read depth (DP) to the INFO, FORMAT and sample columns, in place."""
    depth = sum_read_depth(fields[SAMPLE_COLUMN])
    fields[INFO_COLUMN] += f";DP={depth}"
    fields[FORMAT_COLUMN] += ":DP"
    fields[SAMPLE_COLUMN] += f":{depth}"


def filter_vcf(input_vcf, output_vcf):
    """Write a copy of the VCF with N alleles removed and DP added."""
    dp_header_written = False

    with open(input_vcf) as vcf_in, open(output_vcf, "w") as vcf_out:
        for line_number, line in enumerate(vcf_in, start=1):
            # Header lines are copied as they are. The DP declaration goes in front of
            # the first FORMAT header line, or in front of the column header line if
            # the VCF has no FORMAT header lines.
            if line.startswith("#"):
                is_format_or_column_header = line.startswith(("##FORMAT", "#CHROM"))
                if not dp_header_written and is_format_or_column_header:
                    vcf_out.write(DP_HEADER_LINE)
                    dp_header_written = True
                vcf_out.write(line)
                continue

            if not line.strip():
                continue

            fields = line.strip().split("\t")
            if len(fields) <= SAMPLE_COLUMN:
                raise ValueError(
                    f"{input_vcf}, line {line_number}: expected at least "
                    f"{SAMPLE_COLUMN + 1} tab-separated columns, got {len(fields)}"
                )

            # nothing is left of a record that only has N alleles
            if set(fields[ALT_COLUMN].split(",")) == {"N"}:
                continue

            remove_n_alleles(fields)
            add_read_depth(fields)
            vcf_out.write("\t".join(fields) + "\n")


def main():
    args = parse_args()
    filter_vcf(args.input_vcf, args.output_vcf)


if __name__ == "__main__":
    main()
