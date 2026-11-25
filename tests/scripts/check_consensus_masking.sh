#!/usr/bin/env bash
# Check that every N in a filtered duplex consensus BAM (FilterConsensusReads) is explained by the per-base tags
# of the same record at the same index, and that no kept base meets a masking rule. This only holds when the tags
# are in the orientation of SEQ (the converse check is the one that catches mirrored tags: the depth falls at both
# ends of a read, so a mirrored index often explains an N too). A base is masked when, at its index:
#   - the duplex depth cd (ad + bd when the record has no cd) is below the first --min-reads value, the higher of the single-strand depths ad/bd
#     below the second or the lower one below the third
#   - an error rate (ce/cd, ae/ad, be/bd) is above --max-base-error-rate
#   - the consensus base quality (taken from the unfiltered BAM, as masked bases get quality 2) is below
#     --min-base-quality, or the base was already N there
# Needs samtools. Usage:
#   check_consensus_masking.sh <filtered.bam> <unfiltered.bam> [min_reads="4 2 2"] [min_base_quality=50] [max_base_error_rate=0.1]
# The thresholds are the pipeline defaults (--min_reads, --min_base_quality, --max_base_error_rate).
set -euo pipefail

filtered=$1
unfiltered=$2
min_reads=${3:-4 2 2}
min_baseq=${4:-50}
max_error=${5:-0.1}

# the unfiltered records (bases and qualities), then the filtered ones; records are matched by name and read number
{ samtools view "$unfiltered" | awk -F'\t' '{ print "U\t" $1 "\t" int($2 / 64) % 4 "\t" $10 "\t" $11 }'
    samtools view "$filtered" | awk -F'\t' '{ print "F\t" $0 }'
} | awk -F'\t' -v min_reads="$min_reads" -v min_baseq="$min_baseq" -v max_error="$max_error" '
BEGIN {
    split(min_reads, m, " "); min_duplex = m[1]; min_ab = m[2]; min_ba = m[3]
    for (q = 33; q < 127; q++) ord[sprintf("%c", q)] = q - 33
}
$1 == "U" { key = $2 SUBSEP $3; useq[key] = $4; uqual[key] = $5; next }
{
    name = $2; key = name SUBSEP int($3 / 64) % 4; seq = $11
    if (!(key in useq)) { missing++; next }
    delete tags
    for (i = 13; i <= NF; i++) {
        t = substr($i, 1, 2)
        if (t ~ /^(ad|ae|bd|be|cd|ce)$/) { n = split(substr($i, 8), v, ","); for (j = 1; j <= n; j++) tags[t, j] = v[j] }
    }
    for (j = 1; j <= length(seq); j++) {
        is_n = substr(seq, j, 1) == "N"
        ad = tags["ad", j]; bd = tags["bd", j]
        # duplex depth and errors: cd/ce when present, otherwise (as fgbio does) the sums of the two strands
        cd = (("cd", j) in tags) ? tags["cd", j] : ad + bd; ce = (("ce", j) in tags) ? tags["ce", j] : tags["ae", j] + tags["be", j]
        hi = ad > bd ? ad : bd; lo = ad > bd ? bd : ad
        reason = ""
        if (substr(useq[key], j, 1) == "N") reason = "N before filtering"
        else if (cd < min_duplex || hi < min_ab || lo < min_ba) reason = "depth"
        else if ((cd > 0 && ce / cd > max_error) || (ad > 0 && tags["ae", j] / ad > max_error) || (bd > 0 && tags["be", j] / bd > max_error)) reason = "error rate"
        else if (ord[substr(uqual[key], j, 1)] < min_baseq) reason = "base quality"
        if (!is_n) {
            # and the converse: a base that is kept must not meet any of the rules
            if (reason != "") { kept++; if (kept <= 5) printf "base not masked although %s: %s (flag %s) at %d\n", reason, name, $3, j }
            continue
        }
        total++
        if (reason == "") { unexplained++; if (unexplained <= 5) printf "unexplained N: %s (flag %s) at %d: ad=%s bd=%s cd=%s ae=%s be=%s ce=%s\n", name, $3, j, ad, bd, cd, tags["ae", j], tags["be", j], ce }
        else reasons[reason]++
    }
}
END {
    printf "%d Ns", total
    for (r in reasons) printf ", %d %s", reasons[r], r
    printf ", %d unexplained; %d kept bases that meet a masking rule", unexplained, kept
    if (missing) printf " (%d filtered records not found in the unfiltered BAM)", missing
    print ""
    exit (unexplained > 0 || kept > 0 || missing > 0)
}'
