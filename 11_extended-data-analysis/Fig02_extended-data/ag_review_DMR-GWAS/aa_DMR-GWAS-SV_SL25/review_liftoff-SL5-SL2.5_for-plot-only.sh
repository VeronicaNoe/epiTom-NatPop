#!/usr/bin/env bash
set -euo pipefail

# ----------------------------
# Inputs
# ----------------------------
SV_TSV="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/bb_merge-results/sv.id2coord.tsv"
REF_FA="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/SL5.0.fasta"
TGT_FA="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/ITAG2.4_genomic.fasta"

OUTDIR="/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS"
mkdir -p "$OUTDIR"

REF_GFF="$OUTDIR/sv_markers_SL5.gff3"
LIFTED_GFF="$OUTDIR/sv_markers_SL25_liftoff.gff3"
UNMAPPED="$OUTDIR/sv_markers_SL25_unmapped.txt"
OUT_TSV="$OUTDIR/sv_id2coord_sl25.tsv"
OUT_TSV_PRIMARY="$OUTDIR/sv_id2coord_sl25_primary.tsv"

THREADS=16
HALF_WINDOW=50   # lifts a 101-bp interval centered on the SV position

# ----------------------------
# Build mini GFF3 in SL5 coordinates
# seqids must match SL5 FASTA exactly: 0,1,2,...,12
# ----------------------------
awk -v halfw="$HALF_WINDOW" '
BEGIN {
  FS = OFS = "\t"
  print "##gff-version 3"
}
{
  sv  = $1
  chr = $2
  pos = $3

  if (sv == "" || chr == "" || pos == "") next
  if (pos !~ /^[0-9]+$/) next

  start = pos - halfw
  if (start < 1) start = 1
  end = pos + halfw

  gene_id = sv
  tx_id   = sv ".mrna"
  exon_id = sv ".exon1"

  print chr, "SV", "gene", start, end, ".", "+", ".", "ID=" gene_id ";Name=" gene_id
  print chr, "SV", "mRNA", start, end, ".", "+", ".", "ID=" tx_id ";Parent=" gene_id ";Name=" tx_id
  print chr, "SV", "exon", start, end, ".", "+", ".", "ID=" exon_id ";Parent=" tx_id
}
' "$SV_TSV" > "$REF_GFF"

echo "Built source GFF: $REF_GFF"

# ----------------------------
# Liftoff: SL5 -> SL2.5
# target fasta first, then reference fasta
# ----------------------------
liftoff \
  -g "$REF_GFF" \
  -o "$LIFTED_GFF" \
  -u "$UNMAPPED" \
  -p "$THREADS" \
  "$TGT_FA" \
  "$REF_FA"

echo "Liftoff output: $LIFTED_GFF"
echo "Unmapped list:  $UNMAPPED"

# ----------------------------
# Extract lifted coordinates back to TSV
# Keep the lifted midpoint as the projected marker position
# Output columns:
# sv_id   chr_sl25   pos_sl25   seqid_sl25   start_sl25   end_sl25
# ----------------------------
awk '
BEGIN {
  FS = OFS = "\t"
}
$0 !~ /^#/ && $3 == "gene" {
  seqid = $1
  start = $4
  end   = $5
  attrs = $9

  id = ""
  n = split(attrs, a, ";")
  for (i = 1; i <= n; i++) {
    if (a[i] ~ /^ID=/) {
      sub(/^ID=/, "", a[i])
      id = a[i]
    }
  }

  if (id != "") {
    pos = int((start + end) / 2)

    chr = seqid
    sub(/^SL2\.50ch/, "", chr)
    sub(/^ch/, "", chr)

    # if numeric, zero-pad to 2 digits
    if (chr ~ /^[0-9]+$/) {
      chr = sprintf("%02d", chr + 0)
    }

    print id, chr, pos, seqid, start, end
  }
}
' "$LIFTED_GFF" > "$OUT_TSV"

echo "Lifted TSV: $OUT_TSV"

# ----------------------------
# Optional: keep only primary chr01-chr12
# (drop chr00 or anything else)
# ----------------------------
awk 'BEGIN{FS=OFS="\t"} $2 ~ /^(01|02|03|04|05|06|07|08|09|10|11|12)$/ {print}' \
  "$OUT_TSV" > "$OUT_TSV_PRIMARY"

echo "Primary-only TSV: $OUT_TSV_PRIMARY"

# ----------------------------
# Quick summary
# ----------------------------
echo
echo "Counts:"
echo -n "Input SV markers:      "
wc -l < "$SV_TSV"
echo -n "Lifted TSV rows:       "
wc -l < "$OUT_TSV"
echo -n "Primary-only TSV rows: "
wc -l < "$OUT_TSV_PRIMARY"
echo -n "Unmapped rows:         "
wc -l < "$UNMAPPED"
