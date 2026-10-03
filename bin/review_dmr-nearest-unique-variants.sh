#!/usr/bin/env bash
DMR="$( cat $1 | cut -d'_' -f1)"
CHR="$( cat $1 | cut -d'_' -f2)"
CHRNUM=${CHR#chr}
DMR_DIR="/mnt/disk2/vibanez/otherAnalysis/07_SV-from-graphpangenome/aa_dmrs-SL5-bed"
VAR_DIR="/mnt/disk2/vibanez/otherAnalysis/10_unique-SV-from-graphpangenome/ab_unique-graphpan-SV-to-bed"
OUT_DIR="/mnt/disk2/vibanez/otherAnalysis/10_unique-SV-from-graphpangenome/ac_dmr-nearest-variants"
GENOME="/mnt/disk2/vibanez/otherAnalysis/07_SV-from-graphpangenome/SL5.chrom.sizes"


A="$DMR_DIR/${DMR}_dmrs.SL5.bed"
B="$VAR_DIR/${CHR}.SVLEN50.bed.gz"

# Subset DMRs to this chromosome + sort with the genome order
A_chr_sorted=$(mktemp)
B_sorted=$(mktemp)
trap 'rm -f "$A_chr_sorted" "$B_sorted"' EXIT

awk -v c="$CHRNUM" 'BEGIN{FS=OFS="\t"} $1==c {print}' "$A" \
| bedtools sort -g "$GENOME" -i - > "$A_chr_sorted"

zcat "$B" \
| bedtools sort -g "$GENOME" -i - > "$B_sorted"

bedtools closest -sorted -g "$GENOME" \
  -a "$A_chr_sorted" \
  -b "$B_sorted" \
  -d -t first \
| gzip -c > "$OUT_DIR/${DMR}_${CHR}.closest.tsv.gz"
