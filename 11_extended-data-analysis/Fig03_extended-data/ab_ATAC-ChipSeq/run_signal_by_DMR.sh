#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# Compute ChIP-seq and ATAC-seq bigWig signal over DMRs
# For both C-DMR and CG-DMR
# =========================================================

cd /mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq

# If needed, activate the env that contains deepTools
# conda activate deeptools
# module load deeptools

THREADS=8

DMRSETS=("CG-DMR" "C-DMR")

OUTDIR="aa_prepare-public-data/signal_by_DMR"
mkdir -p "${OUTDIR}"

# ---------------------------------------------------------
# BigWig files
# ---------------------------------------------------------

BW_DIR="aa_prepare-public-data/atlas_bw"

CHIP_BWS=(
  "${BW_DIR}/GSM7844494_coverage_M82_H3K4me3_rep1_S1.bigwig"
  "${BW_DIR}/GSM7844495_coverage_M82_H3K9ac_rep1_S4.bigwig"
  "${BW_DIR}/GSM7844499_coverage_M82_H3K27ac_rep1_S7.bigwig"
  "${BW_DIR}/GSM7844501_coverage_M82_H3K36me3_S6.bigwig"
  "${BW_DIR}/GSM7844511_coverage_M82_H3K9me2_S10.bigwig"
  "${BW_DIR}/GSM7844490_coverage_M82_Pol2_rep1_S10.bigwig"
  "${BW_DIR}/GSM7844510_coverage_M82_H3K27me3_rep1_S1.bigwig"
)

CHIP_LABELS=(
  "H3K4me3"
  "H3K9ac"
  "H3K27ac"
  "H3K36me3"
  "H3K9me2"
  "Pol2"
  "H3K27me3"
)

ATAC_BWS=(
  "${BW_DIR}/GSM8217933_ATAC_M82_rep1_S8_R1.s3norm.bigwig"
  "${BW_DIR}/GSM8217934_ATAC_M82_rep2_S9_R1_001.s3norm.bigwig"
)

ATAC_LABELS=(
  "ATAC_rep1"
  "ATAC_rep2"
)

# ---------------------------------------------------------
# Check files
# ---------------------------------------------------------

for bw in "${CHIP_BWS[@]}" "${ATAC_BWS[@]}"; do
  if [[ ! -s "${bw}" ]]; then
    echo "ERROR: missing bigWig file: ${bw}" >&2
    exit 1
  fi
done

for dmrset in "${DMRSETS[@]}"; do

  BED="aa_prepare-public-data/master/${dmrset}.M82.bed4"

  if [[ ! -s "${BED}" ]]; then
    echo "ERROR: missing BED file: ${BED}" >&2
    exit 1
  fi

  echo "======================================================"
  echo "Processing ${dmrset}"
  echo "BED: ${BED}"
  echo "======================================================"

  # -------------------------------------------------------
  # ChIP signal over DMRs
  # -------------------------------------------------------

  multiBigwigSummary BED-file \
    --BED "${BED}" \
    --bwfiles "${CHIP_BWS[@]}" \
    --labels "${CHIP_LABELS[@]}" \
    --outFileName "${OUTDIR}/${dmrset}.ChIP_signal_by_DMR.npz" \
    --outRawCounts "${OUTDIR}/${dmrset}.ChIP_signal_by_DMR.tsv" \
    --numberOfProcessors "${THREADS}"

  # -------------------------------------------------------
  # ATAC signal over DMRs
  # Optional, but useful because it makes ATAC and ChIP
  # directly comparable: both are bigWig means over DMRs.
  # -------------------------------------------------------

  multiBigwigSummary BED-file \
    --BED "${BED}" \
    --bwfiles "${ATAC_BWS[@]}" \
    --labels "${ATAC_LABELS[@]}" \
    --outFileName "${OUTDIR}/${dmrset}.ATAC_signal_by_DMR.npz" \
    --outRawCounts "${OUTDIR}/${dmrset}.ATAC_signal_by_DMR.tsv" \
    --numberOfProcessors "${THREADS}"

  echo "Finished ${dmrset}"

done

echo "All done."
