#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

TAR="/mnt/disk2/vibanez/otherAnalysis/00_public-data/GSE245529_RAW.tar"
OUT="/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq/aa_prepare-public-data/atlas_bw"
mkdir -p "$OUT"

# Try to extract common M82 bigwigs for these marks (pattern-based)
# You can re-run safely; tar will just overwrite same files.
for pat in "coverage_M82_H3K4me3" "coverage_M82_H3K4me3_rep" \
           "coverage_M82_H3K36me3" "coverage_M82_H3K36me3_rep" \
           "coverage_M82_H3K27ac" "coverage_M82_H3K9ac" "coverage_M82_H3K9me2" \
           "ATAC_M82_rep1" "ATAC_M82_rep2" "coverage_M82_Pol2_rep1" \
	   "coverage_M82_H3K27me3_rep1"
do
  hits=$(tar -tf "$TAR" | grep -F "$pat" | grep -E '\.bigwig$' || true)
  if [[ -n "$hits" ]]; then
    echo "$hits" | while read -r f; do
      tar -xf "$TAR" -C "$OUT" "$f"
    done
  fi
done

ls -lh "$OUT" | head
