#!/usr/bin/env python3
"""
Compare the outputs of the Snakemake pipeline (0.1.0) and the Nextflow pipeline (0.2.0) on the same input.

Usage: compare_parity.py SNAKEMAKE_RESULTS NEXTFLOW_OUTDIR [SAMPLE]

Every comparison prints one line: IDENTICAL, or DIFFERENT with the numbers. Reads are compared field by
field (position, MAPQ, CIGAR, sequence, qualities, mate fields and selected tags). The duplicate flag
(0x400) is masked, because the two pipelines mark duplicates in different BAMs.
"""

import csv
import glob
import gzip
import io
import re
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict

DUP_FLAG = 0x400
KEY_TAGS = ("RX", "ZA", "ZB", "MC", "MD", "NM", "AS", "MQ", "UQ", "RG", "MI", "cD", "cM", "cE")


def samtools_view(bam, header=False):
    cmd = ["samtools", "view"] + (["-H"] if header else []) + [bam]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines()


def read_records(bam):
    """Return a list of dicts, one per alignment record."""
    records = []
    for line in samtools_view(bam):
        f = line.split("\t")
        tags = {}
        for t in f[11:]:
            name, _, rest = t.partition(":")
            tags[name] = rest.partition(":")[2]
        records.append(
            {
                "qname": f[0],
                "flag": int(f[1]) & ~DUP_FLAG,
                "dup": bool(int(f[1]) & DUP_FLAG),
                "rname": f[2],
                "pos": f[3],
                "mapq": f[4],
                "cigar": f[5],
                "rnext": f[6],
                "pnext": f[7],
                "tlen": f[8],
                "seq": f[9],
                "qual": f[10],
                "tags": tags,
            }
        )
    return records


def report(label, ok, detail=""):
    print(f"{'IDENTICAL' if ok else 'DIFFERENT':<10} {label}" + (f"  [{detail}]" if detail else ""))


def compare_reads_by_name(label, a_bam, b_bam):
    """Compare two BAMs record by record, matching records on (read name, flag)."""
    a = defaultdict(list)
    b = defaultdict(list)
    for r in read_records(a_bam):
        a[(r["qname"], r["flag"])].append(r)
    for r in read_records(b_bam):
        b[(r["qname"], r["flag"])].append(r)
    only_a = set(a) - set(b)
    only_b = set(b) - set(a)
    n_a = sum(len(v) for v in a.values())
    n_b = sum(len(v) for v in b.values())
    field_diffs = Counter()
    tag_diffs = Counter()
    examples = {}
    for key in set(a) & set(b):
        for ra, rb in zip(a[key], b[key]):
            for fld in ("rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual"):
                if ra[fld] != rb[fld]:
                    field_diffs[fld] += 1
                    examples.setdefault(fld, (key[0], ra[fld], rb[fld]))
            for t in KEY_TAGS:
                if ra["tags"].get(t) != rb["tags"].get(t):
                    tag_diffs[t] += 1
                    examples.setdefault("tag " + t, (key[0], ra["tags"].get(t), rb["tags"].get(t)))
    ok = not only_a and not only_b and not field_diffs and not tag_diffs and n_a == n_b
    report(label, ok, f"records {n_a} vs {n_b}; only in Snakemake {len(only_a)}, only in Nextflow {len(only_b)}")
    for fld, n in field_diffs.most_common():
        print(f"{'':<11}field {fld}: {n} records differ, e.g. {examples[fld]}")
    for t, n in tag_diffs.most_common():
        print(f"{'':<11}tag {t}: {n} records differ, e.g. {examples['tag ' + t]}")
    return ok


def compare_reads_as_multisets(label, a_bam, b_bam, with_position=True):
    """Compare BAMs whose read names differ (consensus reads): compare multisets of the records."""

    def key(r):
        k = (r["flag"], r["seq"], r["qual"])
        if with_position:
            k = (r["rname"], r["pos"], r["cigar"], r["mapq"], r["rnext"], r["pnext"]) + k
        return k

    ca = Counter(key(r) for r in read_records(a_bam))
    cb = Counter(key(r) for r in read_records(b_bam))
    common = sum((ca & cb).values())
    ok = ca == cb
    report(label, ok, f"records {sum(ca.values())} vs {sum(cb.values())}; identical records {common}")
    return ok


def compare_grouping(label, a_bam, b_bam):
    """Compare the UMI families (MI tag) of two grouped BAMs as partitions of the reads."""

    def partitions(bam):
        by_mi = defaultdict(set)
        by_mi_strand = defaultdict(set)
        for r in read_records(bam):
            mi = r["tags"].get("MI")
            if mi is None:
                continue
            ident = (r["qname"], r["flag"] & 0xC0)
            by_mi[mi.split("/")[0]].add(ident)
            by_mi_strand[mi].add(ident)
        return {frozenset(v) for v in by_mi.values()}, {frozenset(v) for v in by_mi_strand.values()}

    pa_mol, pa_strand = partitions(a_bam)
    pb_mol, pb_strand = partitions(b_bam)
    ok = pa_mol == pb_mol and pa_strand == pb_strand
    report(
        label,
        ok,
        f"molecule families {len(pa_mol)} vs {len(pb_mol)} (same {len(pa_mol & pb_mol)}); "
        f"strand families {len(pa_strand)} vs {len(pb_strand)} (same {len(pa_strand & pb_strand)})",
    )
    return ok


def header_summary(bam):
    lines = samtools_view(bam, header=True)
    hd = [x for x in lines if x.startswith("@HD")]
    rg = [x for x in lines if x.startswith("@RG")]
    sq = [x for x in lines if x.startswith("@SQ")]
    return hd, rg, len(sq)


def read_vcf(path):
    meta, recs = [], []
    for line in open(path):
        if line.startswith("##"):
            meta.append(line.rstrip("\n"))
        elif line.startswith("#"):
            continue
        elif line.strip():
            f = line.rstrip("\n").split("\t")
            recs.append(tuple(f[:5]) + (f[7], f[8], f[9] if len(f) > 9 else ""))
    return meta, recs


def compare_vcf(label, a, b):
    ma, ra = read_vcf(a)
    mb, rb = read_vcf(b)
    ok = ra == rb
    report(label, ok, f"records {len(ra)} vs {len(rb)}; identical records {len(set(ra) & set(rb))}")
    return ok


def read_metrics_csv(path):
    return {r["metric"]: float(r["value"]) for r in csv.DictReader(open(path))}


def compare_metrics(label, a, b, rel_tol=1e-9):
    ma, mb = read_metrics_csv(a), read_metrics_csv(b)
    diffs = {
        m: (ma.get(m), mb.get(m))
        for m in sorted(set(ma) | set(mb))
        if ma.get(m) is None or mb.get(m) is None or abs(ma[m] - mb[m]) > rel_tol * max(1.0, abs(ma[m]))
    }
    report(label, not diffs, f"{len(ma)} vs {len(mb)} metrics, {len(diffs)} differ")
    for m, (x, y) in diffs.items():
        print(f"{'':<11}{m}: Snakemake {x} vs Nextflow {y}")
    return not diffs


def sn_lines(path):
    """The summary numbers (SN lines) of samtools stats, by name, without the command line."""
    out = {}
    for line in open(path):
        if line.startswith("SN\t"):
            f = line.rstrip("\n").split("\t")
            out[f[1].rstrip(":")] = f[2]
    return out


def compare_stats(label, a, b):
    sa, sb = sn_lines(a), sn_lines(b)
    diffs = {k: (sa.get(k), sb.get(k)) for k in sorted(set(sa) | set(sb)) if sa.get(k) != sb.get(k)}
    report(label, not diffs, f"{len(sa)} summary numbers, {len(diffs)} differ")
    for k, (x, y) in diffs.items():
        print(f"{'':<11}{k}: Snakemake {x} vs Nextflow {y}")
    return not diffs


def flagstat_numbers(path):
    out = {}
    for line in open(path):
        m = re.match(r"^(\d+) \+ (\d+) (.*)$", line.strip())
        if m:
            out[re.sub(r"\s*\(.*\)$", "", m.group(3))] = int(m.group(1)) + int(m.group(2))
    return out


def compare_flagstat(label, a, b):
    fa, fb = flagstat_numbers(a), flagstat_numbers(b)
    diffs = {k: (fa.get(k), fb.get(k)) for k in sorted(set(fa) | set(fb)) if fa.get(k) != fb.get(k)}
    report(label, not diffs, f"{len(fa)} lines, {len(diffs)} differ")
    for k, (x, y) in diffs.items():
        print(f"{'':<11}{k}: Snakemake {x} vs Nextflow {y}")
    return not diffs


def picard_metrics(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]
    return dict(zip(rows[0], rows[1]))


def compare_picard(label, a, b):
    ma, mb = picard_metrics(a), picard_metrics(b)
    keep = [k for k in ma if k in mb and k != "LIBRARY"]
    diffs = {k: (ma[k], mb[k]) for k in keep if ma[k] != mb[k]}
    report(label, not diffs, f"{len(keep)} metrics, {len(diffs)} differ")
    for k, (x, y) in diffs.items():
        print(f"{'':<11}{k}: Snakemake {x} vs Nextflow {y}")
    return not diffs


def qualimap_numbers(path):
    out = {}
    for line in open(path):
        m = re.match(r"^\s*(?:>>>>>>> )?([A-Za-z][A-Za-z0-9 /()_.%-]*?) = (.+?)\s*$", line)
        if m:
            out[m.group(1).strip()] = m.group(2)
    # the input file name differs and is not a result
    for k in ("bam file", "outfile", "command line", "date", "bam_file"):
        out.pop(k, None)
    return out


def compare_qualimap(label, a, b):
    qa, qb = qualimap_numbers(a), qualimap_numbers(b)
    diffs = {k: (qa.get(k), qb.get(k)) for k in sorted(set(qa) | set(qb)) if qa.get(k) != qb.get(k)}
    skip = ("bam file", "output", "date", "command", "file")
    diffs = {k: v for k, v in diffs.items() if not any(s in k.lower() for s in skip)}
    report(label, not diffs, f"{len(qa)} values, {len(diffs)} differ")
    for k, (x, y) in diffs.items():
        print(f"{'':<11}{k}: Snakemake {x} vs Nextflow {y}")
    return not diffs


def fastqc_basic(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        name = [n for n in z.namelist() if n.endswith("fastqc_data.txt")][0]
        text = z.read(name).decode()
    out = {}
    block = text.split(">>Basic Statistics")[1].split(">>END_MODULE")[0]
    for line in block.splitlines()[2:]:
        k, _, v = line.partition("\t")
        if k in ("Total Sequences", "Sequences flagged as poor quality", "Sequence length", "%GC"):
            out[k] = v
    return out


def compare_fastqc(label, a_zip, b_zip):
    # FastQC 0.13 no longer writes "Sequences flagged as poor quality", so only the rows of both are compared
    fa, fb = fastqc_basic(a_zip), fastqc_basic(b_zip)
    common = set(fa) & set(fb)
    same = all(fa[k] == fb[k] for k in common)
    report(label, same, str({k: fa[k] for k in sorted(common)}) if same else f"{fa} vs {fb}")
    return same


def fastq_screen_row(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#") and not l.startswith("%")]
    return rows[1] if len(rows) > 1 else []


def compare_screen(label, a, b):
    ra, rb = fastq_screen_row(a), fastq_screen_row(b)
    report(label, ra == rb, f"{ra} vs {rb}" if ra != rb else f"{ra[:6]}")
    return ra == rb


def first(pattern):
    hits = sorted(glob.glob(pattern))
    return hits[0] if hits else None


def main():
    s_dir, n_dir = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/")
    sample = sys.argv[3] if len(sys.argv) > 3 else "NanoMB1Rep1"
    print(f"Snakemake: {s_dir}\nNextflow:  {n_dir}\nsample:    {sample}\n")
    results = {}

    def both(label, s_path, n_path, fn, *args):
        if not s_path or not n_path:
            print(f"{'MISSING':<10} {label}  [Snakemake: {s_path}; Nextflow: {n_path}]")
            results[label] = False
            return
        results[label] = fn(label, s_path, n_path, *args)

    print("== Alignment and UMI-tagged reads")
    s_mapped = f"{s_dir}/mapped_bams/{sample}_mapped_merged.bam"
    n_mapped = f"{n_dir}/mapped_bams/{sample}_mapped_merged.bam"
    sh, nh = header_summary(s_mapped), header_summary(n_mapped)
    print(f"           header: Snakemake @RG {sh[1]} | Nextflow @RG {nh[1]} | @SQ {sh[2]} vs {nh[2]}")
    both("mapped BAM (records)", s_mapped, n_mapped, compare_reads_by_name)

    print("\n== Duplex consensus")
    both("grouped by UMI (UMI families)", f"{s_dir}/grouped_by_umi/{sample}_grouped.bam", f"{n_dir}/grouped_by_umi/{sample}_grouped.bam", compare_grouping)
    both("grouped by UMI (records)", f"{s_dir}/grouped_by_umi/{sample}_grouped.bam", f"{n_dir}/grouped_by_umi/{sample}_grouped.bam", compare_reads_by_name)
    for stem, pos in (
        ("consensus_unmapped", False),
        ("consensus_mapped_merged", True),
        ("consensus_mapped_merged_filtered", True),
        ("consensus_mapped_merged_filtered_clipped", True),
    ):
        both(f"{stem} (reads as multisets)", f"{s_dir}/consensus/{sample}_{stem}.bam", f"{n_dir}/consensus/{sample}_{stem}.bam", compare_reads_as_multisets, pos)

    print("\n== Variants")
    both("nvc VCF", f"{s_dir}/variants/{sample}.vcf", f"{n_dir}/variants/{sample}.vcf", compare_vcf)
    both("filtered VCF", f"{s_dir}/variants/{sample}_filtered.vcf", f"{n_dir}/variants/{sample}_filtered.vcf", compare_vcf)

    print("\n== Duplex metrics")
    both("calc-duplex-metrics CSV", f"{s_dir}/QC/duplex_metrics/{sample}_metrics.csv", f"{n_dir}/QC/duplex_metrics/{sample}_metrics.csv", compare_metrics)

    print("\n== samtools, mapped sample")
    both("samtools stats (summary numbers)", f"{s_dir}/QC/samtools/stats/{sample}.txt", f"{n_dir}/QC/samtools/stats/{sample}.stats", compare_stats)
    both("samtools flagstat", f"{s_dir}/QC/samtools/flagstat/{sample}.txt", f"{n_dir}/QC/samtools/flagstat/{sample}.flagstat", compare_flagstat)
    print("== samtools, consensus")
    both("samtools stats consensus", f"{s_dir}/QC/consensus/samtools/stats/{sample}.txt", f"{n_dir}/QC/consensus/{sample}_consensus.stats", compare_stats)
    both("samtools flagstat consensus", f"{s_dir}/QC/consensus/samtools/flagstat/{sample}.txt", f"{n_dir}/QC/consensus/{sample}_consensus.flagstat", compare_flagstat)

    print("\n== Picard MarkDuplicates")
    both("MarkDuplicates metrics", f"{s_dir}/QC/mark_duplicates/{sample}.txt", first(f"{n_dir}/QC/mark_duplicates/*.metrics.txt"), compare_picard)

    print("\n== Qualimap")
    both("Qualimap mapped", f"{s_dir}/QC/qualimap/{sample}/genome_results.txt", f"{n_dir}/QC/qualimap/{sample}/genome_results.txt", compare_qualimap)
    both("Qualimap consensus", f"{s_dir}/QC/consensus/qualimap/{sample}/genome_results.txt", f"{n_dir}/QC/consensus/{sample}_consensus/genome_results.txt", compare_qualimap)

    print("\n== Raw read QC")
    for r in ("R1", "R2"):
        both(f"FastQC {r} (basic statistics)", first(f"{s_dir}/QC/fastQC/{sample}_{r}_fastqc.zip"), first(f"{n_dir}/QC/fastqc/{sample}_{r[1]}_fastqc.zip"), compare_fastqc)
        both(f"FastQ Screen {r}", first(f"{s_dir}/QC/fastqScreen/{sample}_{r}_screen.txt"), first(f"{n_dir}/QC/fastqscreen/*_{r}_screen.txt"), compare_screen)

    print("\n== MultiQC reports present")
    print(f"           Snakemake main: {bool(first(s_dir + '/QC/multiQC/multiqc_report.html'))}, consensus: {bool(first(s_dir + '/QC/consensus/multiQC/multiqc_report.html'))}")
    print(f"           Nextflow single report: {bool(first(n_dir + '/QC/multiqc/multiqc_report.html'))}")

    n_ok = sum(1 for v in results.values() if v)
    print(f"\nSummary: {n_ok} of {len(results)} comparisons identical")


if __name__ == "__main__":
    main()
