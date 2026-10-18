#!/usr/bin/env bash
set -euo pipefail
SAMPLE="$(cat "$1")"
THREADS="${THREADS:-24}"
DEVICE="cuda:0"

DORADO="/opt/dorado-0.9.6/dorado-0.9.6-linux-x64/bin/dorado"
MODEL="/mnt/ssd123/vibanez/19_ont-mapping/tmp/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
BARCODE_KIT="${BARCODE_KIT:-SQK-NBD114-24}"
POD5_DIR="/mnt/ssd123/vibanez/19_ont-mapping/aa_reads"
REF="/mnt/ssd123/vibanez/19_ont-mapping/ab_run/cervil.polished.fasta"
OUT_DIR="/mnt/ssd123/vibanez/19_ont-mapping/ab_run/${SAMPLE}"

mkdir -p "${OUT_DIR}"

VAR_ALL_BAM="${OUT_DIR}/${SAMPLE}.vars.all.bam"
MOD_ALL_BAM="${OUT_DIR}/${SAMPLE}.mods.all.bam"

VAR_DEMUX_DIR="${OUT_DIR}/vars_demux"
MOD_DEMUX_DIR="${OUT_DIR}/mods_demux"

VAR_SORTED_DIR="${OUT_DIR}/vars_by_barcode"
MOD_SORTED_DIR="${OUT_DIR}/mods_by_barcode"

echo "[1/4] basecall + align + barcode classification for variants/SVs"

#"${DORADO}" basecaller "${MODEL}" "${POD5_DIR}" \
#  --device "${DEVICE}" \
#  --reference "${REF}" \
#  --kit-name "${BARCODE_KIT}" \
#  > "${VAR_ALL_BAM}"

#samtools quickcheck -v "${VAR_ALL_BAM}"
echo "[2/4] demux variants/SVs"
#rm -rf "${VAR_DEMUX_DIR}" "${VAR_SORTED_DIR}"
mkdir -p "${VAR_DEMUX_DIR}" "${VAR_SORTED_DIR}"

"${DORADO}" demux "${VAR_ALL_BAM}" \
	--no-classify \
	--no-trim \
	--output-dir "${VAR_DEMUX_DIR}"

#for f in "${VAR_DEMUX_DIR}"/*.bam; do
#  echo -e "$(basename "$f")\t$(samtools view -c "$f")"
#done | sort -k2,2nr > "${OUT_DIR}/${SAMPLE}.vars.demux.counts.tsv"

for bam in "${VAR_DEMUX_DIR}"/*.bam; do
  base="$(basename "$bam" .bam)"
  if [[ "$base" == *"unclassified"* ]]; then
    echo "Skipping unclassified: ${base}"
    continue
  fi
  out="${VAR_SORTED_DIR}/${base}.vars.sorted.bam"
  samtools sort -@ "${THREADS}" -o "${out}" "${bam}"
  samtools index -@ "${THREADS}" "${out}"
done

echo "[3/4] basecall + align + barcode classification for modified bases"

#"${DORADO}" basecaller sup "${POD5_DIR}" \
#  --device "${DEVICE}" \
#  --reference "${REF}" \
#  --kit-name "${BARCODE_KIT}" \
#  --modified-bases 5mC_5hmC \
#  --modified-bases-threshold 0 \
#  > "${MOD_ALL_BAM}"
#samtools quickcheck -v "${MOD_ALL_BAM}"

echo "[4/4] demux modified-base BAM"
#rm -rf "${MOD_DEMUX_DIR}" "${MOD_SORTED_DIR}"
mkdir -p "${MOD_DEMUX_DIR}" "${MOD_SORTED_DIR}"

"${DORADO}" demux "${MOD_ALL_BAM}" \
        --no-classify \
        --no-trim \
	--output-dir "${MOD_DEMUX_DIR}"

#for f in "${MOD_DEMUX_DIR}"/*.bam; do
#  echo -e "$(basename "$f")\t$(samtools view -c "$f")"
#done | sort -k2,2nr > "${OUT_DIR}/${SAMPLE}.mods.demux.counts.tsv"

for bam in "${MOD_DEMUX_DIR}"/*.bam; do
  base="$(basename "$bam" .bam)"

  if [[ "$base" == *"unclassified"* ]]; then
    echo "Skipping unclassified: ${base}"
    continue
  fi

  out="${MOD_SORTED_DIR}/${base}.mods.sorted.bam"
  samtools sort -@ "${THREADS}" -o "${out}" "${bam}"
  samtools index -@ "${THREADS}" "${out}"
done

echo "[done] ${SAMPLE}"

#echo "Variant demux counts:"
#cat "${OUT_DIR}/${SAMPLE}.vars.demux.counts.tsv"

#echo "Modified-base demux counts:"
#cat "${OUT_DIR}/${SAMPLE}.mods.demux.counts.tsv"
