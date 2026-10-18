#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./barcode_bam_summary.sh /path/to/barcode_bam_dir > barcode_summary.tsv
#
# Example:
#   ./barcode_bam_summary.sh pool-1_test > pool-1_test.summary.tsv
#breadth_pct is the percent of reference bases covered by at least one read.
BAMDIR="${1:-.}"

printf "barcode\tbam\treads\tmapped_reads\ttotal_bases\tmean_read_length_bp\tcovered_bases\tref_bases\tbreadth_pct\tmean_depth\n"

for bam in "$BAMDIR"/*barcode*.bam "$BAMDIR"/*unclassified.bam; do
    [[ -e "$bam" ]] || continue

    base="$(basename "$bam")"

    # Extract barcode name from filename
    barcode="$(echo "$base" | sed -E 's/^.*_(SQK-[^_]+_barcode[0-9]+|unclassified)\.bam$/\1/')"
    if [[ "$barcode" == "$base" ]]; then
        barcode="unknown"
    fi

    # Read counts
    reads="$(samtools view -c "$bam")"
    mapped_reads="$(samtools view -c -F 4 "$bam")"

    # Total bases from read sequence length
    total_bases="$(samtools view "$bam" | awk '{s+=length($10)} END{print s+0}')"

    # Mean read length
    mean_len="$(awk -v b="$total_bases" -v n="$reads" 'BEGIN{if(n>0) printf "%.2f", b/n; else print 0}')"

    # Coverage summary across all reference contigs
    # samtools coverage columns:
    # rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
    read -r covered_bases ref_bases breadth_pct mean_depth < <(
        samtools coverage "$bam" 2>/dev/null | \
        awk '
            BEGIN{covb=0; reflen=0; depthsum=0}
            /^#/ {next}
            NF>=7 {
                len = $3 - $2 + 1
                covb += $5
                reflen += len
                depthsum += $7 * len
            }
            END{
                if(reflen>0){
                    printf "%.0f %.0f %.6f %.6f\n", covb, reflen, 100*covb/reflen, depthsum/reflen
                } else {
                    print "0 0 0 0"
                }
            }'
    )

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$barcode" "$base" "$reads" "$mapped_reads" "$total_bases" "$mean_len" \
        "$covered_bases" "$ref_bases" "$breadth_pct" "$mean_depth"
done | sort -k1,1V
