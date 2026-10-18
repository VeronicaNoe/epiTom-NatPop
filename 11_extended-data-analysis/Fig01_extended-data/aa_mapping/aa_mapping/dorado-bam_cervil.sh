#!/usr/bin/env bash
SAMPLE="$(cat "$1")"
THREADS="${THREADS:-24}"
DEVICE="${DEVICE:-cuda:0}"
DORADO="/opt/dorado-0.9.6/dorado-0.9.6-linux-x64/bin/dorado"

MODEL="/mnt/ssd123/vibanez/19_ont-mapping/tmp/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
POD5_DIR="/mnt/ssd123/vibanez/19_ont-mapping/aa_reads/"
REF="/mnt/ssd123/vibanez/19_ont-mapping/ab_run/cervil.polished.fasta"
OUT_DIR="/mnt/ssd123/vibanez/19_ont-mapping/ab_run/${SAMPLE}"

mkdir -p "${OUT_DIR}"

VAR_BAM="${OUT_DIR}/${SAMPLE}.vars.sorted.bam"
MOD_BAM="${OUT_DIR}/${SAMPLE}.mods.sorted.bam"

echo "[1/2] basecall + align for small variants / SVs: ${SAMPLE}"
"${DORADO}" basecaller "${MODEL}" "${POD5_DIR}" \
  --device "${DEVICE}" \
  --reference "${REF}" \
| samtools sort -@ "${THREADS}" -o "${VAR_BAM}" -

samtools index -@ "${THREADS}" "${VAR_BAM}"

echo "[2/2] basecall + align for modified bases: ${SAMPLE}"
"${DORADO}" basecaller "${MODEL}" "${POD5_DIR}" \
  --device "${DEVICE}" \
  --reference "${REF}" \
  --modified-bases 5mC_5hmC \
  --modified-bases-threshold 0 \
| samtools sort -@ "${THREADS}" -o "${MOD_BAM}" -

samtools index -@ "${THREADS}" "${MOD_BAM}"

echo "[done] ${SAMPLE}"
ls -lh "${VAR_BAM}" "${VAR_BAM}.bai" "${MOD_BAM}" "${MOD_BAM}.bai"
