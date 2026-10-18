#!/usr/bin/env bash
set -euo pipefail

SAMPLE="$(cat "$1")"
THREADS="${THREADS:-16}"
MAPQ="${MAPQ:-20}"

INDIR="/mnt/ssd123/vibanez/19_ont-mapping/ab_run/merged_by_sample/mods_bam"
OUTDIR="/mnt/ssd123/vibanez/19_ont-mapping/ac_modkit"
PILEUP_DIR="${OUTDIR}/pileup"
CX_DIR="${OUTDIR}/cx_report"
LOGDIR="${OUTDIR}/logs"
TMPDIR="${OUTDIR}/tmp_filter"

REF="/mnt/ssd123/vibanez/19_ont-mapping/ab_run/cervil.polished.fasta"
MODKIT="/opt/dist_modkit_v0.4.1_cec0a0b/modkit"
SAMTOOLS="samtools"

BAM="${INDIR}/${SAMPLE}.mods.merged.sorted.bam"

PILEUP_BED="${PILEUP_DIR}/${SAMPLE}.pileup_all.bed"
CX_REPORT="${CX_DIR}/${SAMPLE}.CX_report.txt.gz"
LOG="${LOGDIR}/${SAMPLE}.modkit_pileup.log"

mkdir -p "${PILEUP_DIR}" "${CX_DIR}" "${LOGDIR}" "${TMPDIR}"

[[ -s "${BAM}" ]] || { echo "[ERROR] Missing BAM: ${BAM}" >&2; exit 1; }
[[ -s "${BAM}.bai" ]] || { echo "[ERROR] Missing BAM index: ${BAM}.bai" >&2; exit 1; }
[[ -s "${REF}" ]] || { echo "[ERROR] Missing reference: ${REF}" >&2; exit 1; }

TMP_BAM="${TMPDIR}/${SAMPLE}.primary.mapq${MAPQ}.bam"

cleanup() {
  rm -f "${TMP_BAM}" "${TMP_BAM}.bai"
}
trap cleanup EXIT

echo "[filter] ${SAMPLE}: keep primary alignments with MAPQ >= ${MAPQ}" >&2

TOTAL=$("${SAMTOOLS}" view -c "${BAM}")
KEPT=$("${SAMTOOLS}" view -c -q "${MAPQ}" -F 2308 "${BAM}")

echo "[filter] total alignments: ${TOTAL}" >&2
echo "[filter] kept alignments:  ${KEPT}" >&2

"${SAMTOOLS}" view -bh -q "${MAPQ}" -F 2308 "${BAM}" > "${TMP_BAM}"
"${SAMTOOLS}" index -@ "${THREADS}" "${TMP_BAM}"

echo "[modkit] ${SAMPLE}" >&2

rm -f "${PILEUP_BED}" "${CX_REPORT}" "${LOG}"

"${MODKIT}" pileup \
  "${TMP_BAM}" \
  "${PILEUP_BED}" \
  --ref "${REF}" \
  --cpg \
  --motif CHG 0 \
  --motif CHH 0 \
  --threads "${THREADS}" \
  --filter-threshold C:0.70 \
  --edge-filter 5 \
  --log-filepath "${LOG}"

echo "[convert] ${SAMPLE} pileup -> CX report" >&2

awk -F'\t' 'BEGIN{OFS="\t"}
  $4 ~ /^m,/ {
    split($4,a,",");
    ctx=a[2];
    print $1, $2+1, $6, $12, $13, ctx, ctx
  }
' "${PILEUP_BED}" | gzip -c > "${CX_REPORT}"

echo "[done] ${CX_REPORT}" >&2
