#!/usr/bin/env bash
set -euo pipefail

REF="/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/02_polish/cervil.polished.fasta"
PROJ="/work/project/ijpb/vibanez/14_cervil-project/04_gen-var"
DV="${PROJ}/aa_run-deepvariant"
OUT="${PROJ}/ac_compare"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

module load bioinfo/Bcftools/1.21
module load bioinfo/tabix/1.18
module load bioinfo/samtools/1.21
module load devel/Miniforge/Miniforge3


mkdir -p "${OUT}/gvcf_reheader" "${OUT}/logs"

[[ -f "${REF}.fai" ]] || samtools faidx "${REF}"

JOINT_BCF="${OUT}/all-samples.deepvariant.glnexus.bcf"
JOINT_VCF="${OUT}/all-samples.deepvariant.glnexus.vcf.gz"
MASKED="${OUT}/all-samples.deepvariant.glnexus.GQ20.DP10.vcf.gz"
TAGS="${OUT}/all-samples.deepvariant.glnexus.GQ20.DP10.tags.vcf.gz"

rm -f "${OUT}/gvcf_reheader/"*.g.vcf.gz*
rm -f "${OUT}/all-samples.deepvariant.glnexus"*

echo "[reheader gVCFs]"

for S in G0-C3 G0-C6 G0-C9 G5-C1 G5-C2 G5-C3; do
    IN="${DV}/${S}/${S}.deepvariant.g.vcf.gz"
    RH="${OUT}/gvcf_reheader/${S}.deepvariant.g.vcf.gz"
    SAMPLE_TXT="${OUT}/gvcf_reheader/${S}.sample.txt"

    echo "[sample] ${S}"

    [[ -f "${IN}" ]] || { echo "Missing input: ${IN}"; exit 1; }

    echo "${S}" > "${SAMPLE_TXT}"

    bcftools reheader \
        -s "${SAMPLE_TXT}" \
        -o "${RH}" \
        "${IN}"

    tabix -f -p vcf "${RH}"
done

GVCFS=(
  "${OUT}/gvcf_reheader/G0-C3.deepvariant.g.vcf.gz"
  "${OUT}/gvcf_reheader/G0-C6.deepvariant.g.vcf.gz"
  "${OUT}/gvcf_reheader/G0-C9.deepvariant.g.vcf.gz"
  "${OUT}/gvcf_reheader/G5-C1.deepvariant.g.vcf.gz"
  "${OUT}/gvcf_reheader/G5-C2.deepvariant.g.vcf.gz"
  "${OUT}/gvcf_reheader/G5-C3.deepvariant.g.vcf.gz"
)

echo "[joint genotype with GLnexus]"

glnexus_cli \
    --config DeepVariant \
    --threads "${THREADS}" \
    "${GVCFS[@]}" \
    > "${JOINT_BCF}"

bcftools index -f "${JOINT_BCF}"

echo "[convert BCF to VCF.GZ]"

bcftools view \
    -Oz \
    -o "${JOINT_VCF}" \
    "${JOINT_BCF}"

tabix -f -p vcf "${JOINT_VCF}"

echo "[mask genotypes: GQ < 20 or DP < 10]"

bcftools +setGT "${JOINT_VCF}" \
    -Oz \
    -o "${MASKED}" \
    -- \
    -t q \
    -n . \
    -i 'FMT/GQ<20 || FMT/DP<10'

tabix -f -p vcf "${MASKED}"

echo "[fill tags]"

bcftools +fill-tags "${MASKED}" \
    -Oz \
    -o "${TAGS}" \
    -- -t F_MISSING,NS,AN,AC,AF

tabix -f -p vcf "${TAGS}"

echo "[done]"
echo "Joint VCF:  ${JOINT_VCF}"
echo "Masked VCF: ${MASKED}"
echo "Tagged VCF: ${TAGS}"
