#!/usr/bin/env python3
"""Check the sample name cleaning of assets/multiqc_config.yml (extra_fn_clean_exts, regex entries) on the
names the pipeline's files give: the "<sample>_mapped" / "<sample>_consensus" / "<sample>_reads_N" rows have to
come out right, also for sample names that contain parts of the suffixes.

Usage: python3 tests/scripts/check_multiqc_names.py
"""
import re
import sys
from pathlib import Path

config = (Path(__file__).resolve().parents[2] / "assets" / "multiqc_config.yml").read_text()
block = config.split("extra_fn_clean_exts:")[1]
patterns = [p.replace("\\\\", "\\") for p in re.findall(r'type: regex\s+pattern: "([^"]+)"', block)]  # YAML double-quoted escapes
truncates = re.findall(r'type: truncate\s+pattern: "([^"]+)"', block)


def clean(name):
    for t in truncates:
        name = name.split(t)[0]
    for p in patterns:
        name = re.sub(p, "", name)
    return name


cases = {
    "S_mapped_merged_markdup": "S_mapped",
    "S_mapped_merged": "S_mapped",
    # Picard names the sample after its --INPUT, with the extension still on when the cleaning runs
    "S_mapped_merged.bam": "S_mapped",
    "S_mapped_merged_markdup.bam": "S_mapped",
    "S_consensus_mapped_merged_filtered.bam": "S_consensus",
    "S_mapped": "S_mapped",
    "S_consensus": "S_consensus",
    "S_consensus.mosdepth.summary": "S_consensus",
    "S_consensus_mapped_merged_filtered": "S_consensus",
    "S_reads_1": "S_reads_1",
    # sample names that contain parts of the suffixes
    "tumour_merged_A_reads_1": "tumour_merged_A_reads_1",
    "tumour_merged_A_mapped_merged_markdup": "tumour_merged_A_mapped",
    "tumour_merged_A_consensus_mapped_merged_filtered": "tumour_merged_A_consensus",
    "tumour_merged_B_mapped": "tumour_merged_B_mapped",
}

failed = 0
for raw, expected in cases.items():
    got = clean(raw)
    if got != expected:
        failed += 1
        print(f"FAIL {raw}: got {got}, expected {expected}")
print(f"{len(cases) - failed}/{len(cases)} names cleaned as expected")
sys.exit(1 if failed else 0)
