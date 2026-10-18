#!/usr/bin/env bash
set -euo pipefail

REF="/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/02_polish/cervil.polished.fasta"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
OUTDIR_SV="/work/project/ijpb/vibanez/14_cervil-project/04_gen-var/ab_run-sniffles"

module load bioinfo/Bcftools/1.21
module load bioinfo/tabix/1.18

JOINT="${OUTDIR_SV}/all-samples_joint.vcf.gz"

FILTERED="${OUTDIR_SV}/all-samples_joint.PASS.SUP5.vcf.gz"
MASKED="${OUTDIR_SV}/all-samples_joint.PASS.SUP5.maskedGT.vcf.gz"
TAGS="${OUTDIR_SV}/all-samples_joint.PASS.SUP5.maskedGT.tags.vcf.gz"
FINAL="${OUTDIR_SV}/all-samples_joint.PASS.SUP5.maskedGT.max1missing.vcf.gz"
FILTERED="${OUTDIR_SV}/all-samples_joint.PASS.SUP5.vcf.gz"
FILTERED_TAGS="${OUTDIR_SV}/all-samples_joint.PASS.SUP5.tags.vcf.gz"


# -----------------------------
# 1. Joint SV calling
# -----------------------------
sniffles \
  --input \
    "${OUTDIR_SV}/G0-C3.snf" \
    "${OUTDIR_SV}/G0-C6.snf" \
    "${OUTDIR_SV}/G0-C9.snf" \
    "${OUTDIR_SV}/G5-C1.snf" \
    "${OUTDIR_SV}/G5-C2.snf" \
    "${OUTDIR_SV}/G5-C3.snf" \
  --vcf "${JOINT}" \
  --reference "${REF}" \
  --threads "${THREADS}"

tabix -f -p vcf "${JOINT}"

# -----------------------------
# 2. Site-level filter
# Keep confident SV records
# -----------------------------
bcftools view \
  -i 'FILTER="PASS" && INFO/SUPPORT>=5' \
  -Oz -o "${FILTERED}" \
  "${JOINT}"

tabix -f -p vcf "${FILTERED}"
# Add tags to support-filtered file
# This is the file we will use in R

bcftools +fill-tags "${FILTERED}" \
  -Oz -o "${FILTERED_TAGS}" \
  -- -t NS,F_MISSING,AC,AN

tabix -f -p vcf "${FILTERED_TAGS}"

# -----------------------------
# 3. Genotype-level masking
# Weak genotypes become ./.
# -----------------------------
#bcftools +setGT "${FILTERED}" \
#  -- -t q -n . \
#  -i 'FMT/GQ<10 || (FMT/DR+FMT/DV)<5 || (GT!="0/0" && FMT/DV<3)' \
#  | bcftools view -Oz -o "${MASKED}"

#tabix -f -p vcf "${MASKED}"

# -----------------------------
# 4. Add missingness / allele-count tags
# -----------------------------
#bcftools +fill-tags "${MASKED}" \
#  -Oz -o "${TAGS}" \
#  -- -t NS,F_MISSING,AC,AN

#tabix -f -p vcf "${TAGS}"

# -----------------------------
# 5. Variant-level missingness filter
# Keep variants with at most 1 missing sample
# 6 samples total => NS >= 5
# Also remove variants with no ALT allele after masking
# -----------------------------
#bcftools view \
#  -i 'NS>=5 && AC>0' \
#  -Oz -o "${FINAL}" \
#  "${TAGS}"

#tabix -f -p vcf "${FINAL}"

echo "Done."
