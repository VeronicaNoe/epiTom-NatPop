#!/usr/bin/env bash
ROOT="/mnt/disk2/vibanez/otherAnalysis/10_unique-SV-from-graphpangenome"
BEDDIR="$ROOT/ab_unique-graphpan-SV-to-bed"
DMRDIR="/mnt/disk2/vibanez/otherAnalysis/07_SV-from-graphpangenome/aa_dmrs-SL5-bed"
GENOME="$ROOT/SL5.chrom.sizes"
OUT_DIR="$ROOT/ac_intersects-closest"

SVLEN="${1:-50}"

echo "[STEP 1] Build complete all-variant BED track"
zcat "$BEDDIR"/chr*.SVLEN${SVLEN}.bed.gz \
| bedtools sort -g "$GENOME" -i - \
| bgzip -c > "$BEDDIR/allVariants.SVLEN${SVLEN}.bed.gz"
tabix -p bed "$BEDDIR/allVariants.SVLEN${SVLEN}.bed.gz"

for DMR in C-DMR CG-DMR; do
  echo "[STEP 2] Process $DMR"
  A="$DMRDIR/${DMR}_dmrs.SL5.bed"
  A_SORT="$DMRDIR/${DMR}_dmrs.SL5.sorted.bed.gz"
  B="$BEDDIR/allVariants.SVLEN${SVLEN}.bed.gz"
  # sort DMR bed once
  bedtools sort -g "$GENOME" -i "$A" | bgzip -c > "$A_SORT"
  tabix -p bed "$A_SORT"
  # 2A) TRUE overlap with ALL variants
  bedtools intersect -sorted -g "$GENOME" \
    -a <(zcat "$A_SORT") \
    -b <(zcat "$B") \
    -wo \
  | bgzip -c > "$OUT_DIR/${DMR}_intersect-allVariants.tsv.gz"

  # 2B) DMRs with no overlap at all
  bedtools intersect -sorted -g "$GENOME" \
    -a <(zcat "$A_SORT") \
    -b <(zcat "$B") \
    -v \
  | bgzip -c > "$OUT_DIR/${DMR}_dmrs.SL5.noOverlap.bed.gz"
  tabix -p bed "$OUT_DIR/${DMR}_dmrs.SL5.noOverlap.bed.gz"
  # 3) DMRs without overlap
  echo "[STEP 3] Closest for no-overlap $DMR"
  A="$OUT_DIR/${DMR}_dmrs.SL5.noOverlap.bed.gz"
  bedtools closest -sorted -g "$GENOME" \
    -a <(zcat "$A") \
    -b <(zcat "$B") \
    -d -t first \
  | bgzip -c > "$OUT_DIR/${DMR}_closest-allVariantsWOoverlap.tsv.gz"
done
