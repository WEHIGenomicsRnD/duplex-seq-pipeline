#!/usr/bin/env bash
# Check that the per-base consensus tags of fgbio are in the orientation of SEQ on every record of a
# consensus BAM (the realigned, filtered or clipped one):
#   - every per-base tag is as long as SEQ (ClipBam --auto-clip-attributes)
#   - ac (the consensus bases) equals SEQ at every position where SEQ is not N, on both strands
# Needs samtools. Usage: check_consensus_tag_orientation.sh <consensus.bam>
set -euo pipefail

samtools view "$1" | awk -F'\t' '
{
    seq = $10; len = length(seq); rev = int($2 / 16) % 2; ac = ""; len_ok = 1
    for (i = 12; i <= NF; i++) {
        tag = substr($i, 1, 2); value = substr($i, 6)
        if (tag ~ /^(ad|ae|aq|ac|bd|be|bq|bc|cd|ce)$/) {
            # B arrays (ad, ae, ...) are "type,v1,v2,..."; Z strings (ac, aq, ...) are one character per base
            n = substr($i, 4, 1) == "B" ? split(value, parts, ",") - 1 : length(value)
            if (n != len && len_ok) { len_ok = 0; bad_len[rev]++; if (bad_len[rev] == 1) print "tag length differs from SEQ: " $1 " " tag " (" n " vs " len ")" }
        }
        if (tag == "ac") ac = value
    }
    if (ac == "") next
    checked[rev]++
    for (j = 1; j <= len; j++) {
        b = substr(seq, j, 1)
        if (b != "N" && b != substr(ac, j, 1)) { bad_seq[rev]++; if (bad_seq[rev] == 1) print "ac differs from SEQ: " $1 " at " j " (strand " (rev ? "-" : "+") ")"; break }
    }
}
END {
    for (r = 0; r <= 1; r++)
        printf "%s strand: %d records with ac, %d where ac differs from SEQ, %d with a tag of another length than SEQ\n", (r ? "-" : "+"), checked[r], bad_seq[r], bad_len[r]
    exit (bad_seq[0] + bad_seq[1] + bad_len[0] + bad_len[1] > 0)
}'
