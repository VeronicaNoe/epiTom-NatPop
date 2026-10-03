#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

OUTROOT="/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq/ac_annotate-DMRs-M82"
LIFTROOT="/mnt/disk2/vibanez/otherAnalysis/04_liftoff-DMRs-SL2.5-M82"

REG="${OUTROOT}/regions_M82"
G="${REG}/M82.chrom.sizes.chr"
PROM="${REG}/promoter3500.merged.bed3"
GENE="${REG}/genes.merged.bed3"
TE="${REG}/TE.merged.bed3"
INTER="${REG}/intergenic.bed3"

DMRSETS=(C-DMR)   # add CHG-DMR/CHH-DMR if you have them
SEED=1

mkdir -p "${OUTROOT}/dmrs_by_region" "${OUTROOT}/shuffle_by_region"

for DMRSET in "${DMRSETS[@]}"; do
  IN_GZ="${LIFTROOT}/${DMRSET}.M82.hiconf.mapq30.chr.bed4.gz"
  [[ -s "$IN_GZ" ]] || { echo "WARNING: missing $IN_GZ, skipping $DMRSET"; continue; }
  # base BED4
  BASE="${OUTROOT}/dmrs_by_region/${DMRSET}.bed4"
  zcat "$IN_GZ" > "$BASE"
  # exclusive assignment (priority): promoter > gene > TE > intergenic
  # promoter
  bedtools intersect -u -a "$BASE" -b "$PROM" > "${OUTROOT}/dmrs_by_region/${DMRSET}.promoter.bed4"
  # gene (excluding promoter)
  bedtools intersect -u -a "$BASE" -b "$GENE" \
  | bedtools intersect -v -a - -b "$PROM" \
  > "${OUTROOT}/dmrs_by_region/${DMRSET}.gene.bed4"
  # TE (excluding promoter+gene)
  bedtools intersect -u -a "$BASE" -b "$TE" \
  | bedtools intersect -v -a - -b "$PROM" \
  | bedtools intersect -v -a - -b "$GENE" \
  > "${OUTROOT}/dmrs_by_region/${DMRSET}.TE.bed4"
  # intergenic (not overlapping union)
  bedtools intersect -u -a "$BASE" -b "$INTER" \
  > "${OUTROOT}/dmrs_by_region/${DMRSET}.intergenic.bed4"
  # gene_TE and gene_noTE (DMRs in genes that overlap TE vs not)
  bedtools intersect -u -a "${OUTROOT}/dmrs_by_region/${DMRSET}.gene.bed4" -b "$TE" \
  > "${OUTROOT}/dmrs_by_region/${DMRSET}.gene_TE.bed4"
  bedtools intersect -v -a "${OUTROOT}/dmrs_by_region/${DMRSET}.gene.bed4" -b "$TE" \
  > "${OUTROOT}/dmrs_by_region/${DMRSET}.gene_noTE.bed4"
  # Shuffles: within the matching universe bed (same count, length preserved)
  for R in promoter gene TE intergenic gene_TE gene_noTE; do
    IN="${OUTROOT}/dmrs_by_region/${DMRSET}.${R}.bed4"
    [[ -s "$IN" ]] || continue
    case "$R" in
      promoter) U="$PROM" ;;
      gene) U="$GENE" ;;
      TE) U="$TE" ;;
      intergenic) U="$INTER" ;;
      gene_TE) U="${REG}/genes_TE.bed3" ;;
      gene_noTE) U="${REG}/genes_noTE.bed3" ;;
    esac
    OUT="${OUTROOT}/shuffle_by_region/${DMRSET}.${R}.shuffle.bed3"
    bedtools shuffle -i <(cut -f1-3 "$IN") -g "$G" -incl "$U" -chrom -seed "$SEED" \
      > "$OUT"
  done
  echo "== $DMRSET counts =="
  for R in promoter gene TE intergenic gene_TE gene_noTE; do
    f="${OUTROOT}/dmrs_by_region/${DMRSET}.${R}.bed4"
    [[ -s "$f" ]] && echo "$R:" $(wc -l < "$f")
  done
done

echo "Wrote:"
echo "  ${OUTROOT}/dmrs_by_region/"
echo "  ${OUTROOT}/shuffle_by_region/"
