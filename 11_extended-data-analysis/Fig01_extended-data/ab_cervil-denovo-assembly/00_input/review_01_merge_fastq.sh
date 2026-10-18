#!/usr/bin/env bash
PROJECT="/work/project/ijpb/vibanez/14_cervil-project"
NEWPROJ="${PROJECT}/07_cervil-denovo-assembly"
OLD_FASTQ_DIR="${PROJECT}/02_prep/trimmed_fastq"
SAMPLE_LIST="${NEWPROJ}/00_input/sample_list.txt"
OUT_FASTQ="${NEWPROJ}/00_input/cervil_pooled.trimmed.fastq.gz"
TMP_LIST="${NEWPROJ}/00_input/cervil_pooled.inputs.txt"

mkdir -p "${NEWPROJ}/00_input"

: > "${TMP_LIST}"

while read -r SAMPLE; do
    [[ -z "${SAMPLE}" ]] && continue
    FQ="${OLD_FASTQ_DIR}/${SAMPLE}.trimmed.fastq.gz"
    if [[ ! -s "${FQ}" ]]; then
        echo "[ERROR] Missing file: ${FQ}" >&2
        exit 1
    fi
    echo "${FQ}" >> "${TMP_LIST}"
done < "${SAMPLE_LIST}"

echo "[merge] building ${OUT_FASTQ}" >&2
xargs cat < "${TMP_LIST}" > "${OUT_FASTQ}"

echo "[merge] done" >&2
ls -lh "${OUT_FASTQ}" >&2
