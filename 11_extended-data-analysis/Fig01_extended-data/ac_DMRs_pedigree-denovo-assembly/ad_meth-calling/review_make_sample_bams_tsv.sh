#!/usr/bin/env bash
set -euo pipefail

ROOT="/work/project/ijpb/vibanez/cervil_ont/ab_mapping"
OUT="/work/project/ijpb/vibanez/cervil_ont/ac_methylation-calling/sample_bams.tsv"

mkdir -p "$(dirname "$OUT")"

{
  printf "sample\tbams\n"

  printf "G0-C3\t"
  find "${ROOT}"/pool-* -maxdepth 1 -type f -name '*barcode01.bam' | sort | paste -sd, -
  printf ",%s\n" "${ROOT}/G0-C3_s1/calls_2026-03-27_T07-23-43.bam"

  printf "G0-C6\t"
  find "${ROOT}"/pool-* -maxdepth 1 -type f -name '*barcode02.bam' | sort | paste -sd, -
  printf "\n"

  printf "G0-C9\t"
  find "${ROOT}"/pool-* -maxdepth 1 -type f -name '*barcode03.bam' | sort | paste -sd, -
  printf "\n"

  printf "G5-C1\t"
  find "${ROOT}"/pool-* -maxdepth 1 -type f -name '*barcode04.bam' | sort | paste -sd, -
  printf "\n"

  printf "G5-C2\t"
  find "${ROOT}"/pool-* -maxdepth 1 -type f -name '*barcode05.bam' | sort | paste -sd, -
  printf ",%s\n" "${ROOT}/G5-C2_s5/calls_2026-03-26_T13-34-24.bam"

  printf "G5-C3\t"
  find "${ROOT}"/pool-* -maxdepth 1 -type f -name '*barcode06.bam' | sort | paste -sd, -
  printf "\n"
} > "${OUT}"

echo "Wrote: ${OUT}"
