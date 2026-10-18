#!/usr/bin/env bash
set -euo pipefail

OUTROOT="/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/ag_ONT-DMR-genetic-vars-comparison"
BEDDIR="${OUTROOT}/01_bed"
OVDIR="${OUTROOT}/02_overlap"
CLDIR="${OUTROOT}/03_closest"

mkdir -p "${OVDIR}" "${CLDIR}"

command -v bedtools >/dev/null 2>&1 || { echo "bedtools not found in PATH"; exit 1; }

# sort inputs
sort -k1,1 -k2,2n "${BEDDIR}/dmrs.bed" > "${BEDDIR}/dmrs.sorted.bed"

for X in smallvars.all smallvars.informative sv.all sv.informative; do
    sort -k1,1 -k2,2n "${BEDDIR}/${X}.bed" > "${BEDDIR}/${X}.sorted.bed"
done

# overlaps
bedtools intersect -wa -wb \
  -a "${BEDDIR}/dmrs.sorted.bed" \
  -b "${BEDDIR}/smallvars.all.sorted.bed" \
  > "${OVDIR}/dmr_vs_smallvars.all.overlap.tsv"

bedtools intersect -wa -wb \
  -a "${BEDDIR}/dmrs.sorted.bed" \
  -b "${BEDDIR}/smallvars.informative.sorted.bed" \
  > "${OVDIR}/dmr_vs_smallvars.informative.overlap.tsv"

bedtools intersect -wa -wb \
  -a "${BEDDIR}/dmrs.sorted.bed" \
  -b "${BEDDIR}/sv.all.sorted.bed" \
  > "${OVDIR}/dmr_vs_sv.all.overlap.tsv"

bedtools intersect -wa -wb \
  -a "${BEDDIR}/dmrs.sorted.bed" \
  -b "${BEDDIR}/sv.informative.sorted.bed" \
  > "${OVDIR}/dmr_vs_sv.informative.overlap.tsv"

# nearest
bedtools closest -d -t first \
  -a "${BEDDIR}/dmrs.sorted.bed" \
  -b "${BEDDIR}/smallvars.all.sorted.bed" \
  > "${CLDIR}/dmr_vs_smallvars.all.closest.tsv"

bedtools closest -d -t first \
  -a "${BEDDIR}/dmrs.sorted.bed" \
  -b "${BEDDIR}/smallvars.informative.sorted.bed" \
  > "${CLDIR}/dmr_vs_smallvars.informative.closest.tsv"

bedtools closest -d -t first \
  -a "${BEDDIR}/dmrs.sorted.bed" \
  -b "${BEDDIR}/sv.all.sorted.bed" \
  > "${CLDIR}/dmr_vs_sv.all.closest.tsv"

bedtools closest -d -t first \
  -a "${BEDDIR}/dmrs.sorted.bed" \
  -b "${BEDDIR}/sv.informative.sorted.bed" \
  > "${CLDIR}/dmr_vs_sv.informative.closest.tsv"

echo "Done."
