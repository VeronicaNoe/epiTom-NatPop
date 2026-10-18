#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./pairwise_diff.sh SNP
#   ./pairwise_diff.sh DMR BIseq
#   ./pairwise_diff.sh DMR ONT

MARKER="${1:-}"
TECH="${2:-}"

VCF_DIR="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ak_genetic-variants"
DMR_DIR="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/al_DMR-gen-var_denovo-pedigree_comparison"
OUTDIR="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/am_pairwise-divergences"

SAMPLE_LIST="${OUTDIR}/sample_names.txt"
SAMPLES="$(tr '\n' '\t' < "${SAMPLE_LIST}")"
THREADS=80

mkdir -p "${OUTDIR}"

if [[ ! -s "${SAMPLE_LIST}" ]]; then
  echo "ERROR: sample list not found or empty: ${SAMPLE_LIST}" >&2
  exit 1
fi

echo "Samples used by --sample-diff:"
cat "${SAMPLE_LIST}"
echo

# =========================
# SNP pairwise differences
# =========================
if [[ "${MARKER}" == "SNP" ]]; then

  echo "1. SNP pairwise differences"

  SNP_PREFIX="${OUTDIR}/general.SNPs-denovo.GQ20_DP10"

  # remove old output to avoid using stale files
  rm -f "${SNP_PREFIX}".*

  plink2 \
    --vcf "${VCF_DIR}/cervil_SNPs.GQ20.DP10.vcf.gz" \
    --make-bed \
    --vcf-half-call missing \
    --max-alleles 2 \
    --out "${SNP_PREFIX}"

  plink2 \
    --bfile "${SNP_PREFIX}" \
    --double-id \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --sample-diff pairwise ids=${SAMPLES} \
    --out "${SNP_PREFIX}"

  echo "Done SNP."
  exit 0
fi

# =========================
# DMR pairwise differences
# =========================
if [[ "${MARKER}" == "DMR" ]]; then

  if [[ "${TECH}" != "BIseq" && "${TECH}" != "ONT" ]]; then
    echo "ERROR: for DMR, second argument must be BIseq or ONT" >&2
    echo "Usage:" >&2
    echo "  ./pairwise_diff.sh DMR BIseq" >&2
    echo "  ./pairwise_diff.sh DMR ONT" >&2
    exit 1
  fi

  echo "2. ${TECH} C-DMR pairwise differences"

  C_PREFIX="${OUTDIR}/general.${TECH}.C-DMR-denovo"
  rm -f "${C_PREFIX}".*

  plink2 \
    --vcf "${DMR_DIR}/02.4_dmr_binary_states.${TECH}.C-DMR.vcf.gz" \
    --allow-extra-chr \
    --make-pgen \
    --sort-vars \
    --out "${C_PREFIX}"

  echo "Checking sample IDs in ${C_PREFIX}.psam:"
  cat "${C_PREFIX}.psam"
  echo

  plink2 \
    --pfile "${C_PREFIX}" \
    --double-id \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --sample-diff pairwise ids=${SAMPLES} \
    --out "${C_PREFIX}"

  echo "3. ${TECH} CG-DMR pairwise differences"

  CG_PREFIX="${OUTDIR}/general.${TECH}.CG-DMR-denovo"
  rm -f "${CG_PREFIX}".*

  plink2 \
    --vcf "${DMR_DIR}/02.4_dmr_binary_states.${TECH}.CG-DMR.vcf.gz" \
    --allow-extra-chr \
    --make-pgen \
    --sort-vars \
    --out "${CG_PREFIX}"

  echo "Checking sample IDs in ${CG_PREFIX}.psam:"
  cat "${CG_PREFIX}.psam"
  echo

  plink2 \
    --pfile "${CG_PREFIX}" \
    --double-id \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --sample-diff pairwise ids=${SAMPLES} \
    --out "${CG_PREFIX}"

  echo "Done DMR ${TECH}."
  exit 0
fi

# =========================
# Bad input
# =========================
echo "ERROR: first argument must be SNP or DMR" >&2
echo "Usage:" >&2
echo "  ./pairwise_diff.sh SNP" >&2
echo "  ./pairwise_diff.sh DMR BIseq" >&2
echo "  ./pairwise_diff.sh DMR ONT" >&2
exit 1
