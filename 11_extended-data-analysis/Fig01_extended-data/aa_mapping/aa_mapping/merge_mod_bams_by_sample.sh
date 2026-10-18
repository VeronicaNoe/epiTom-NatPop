#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-16}"

ROOT="/mnt/ssd123/vibanez/19_ont-mapping/ab_run"
OUTDIR="${ROOT}/merged_by_sample/mods_bam"
LISTDIR="${ROOT}/merged_by_sample/filelists"

mkdir -p "${OUTDIR}" "${LISTDIR}"

declare -A SAMPLE_TO_BARCODE=(
  [G0-C3]="barcode01"
  [G0-C6]="barcode02"
  [G0-C9]="barcode03"
  [G5-C1]="barcode04"
  [G5-C2]="barcode05"
  [G5-C3]="barcode06"
)

for SAMPLE in G0-C3 G0-C6 G0-C9 G5-C1 G5-C2 G5-C3; do
    BARCODE="${SAMPLE_TO_BARCODE[$SAMPLE]}"
    LIST="${LISTDIR}/${SAMPLE}.mods.bam.list"
    OUTBAM="${OUTDIR}/${SAMPLE}.mods.merged.sorted.bam"

    echo "[sample] ${SAMPLE} ${BARCODE}" >&2

    : > "${LIST}"

    # pooled libraries
    find "${ROOT}/all/mods_by_barcode" \
         "${ROOT}/pool-1/mods_by_barcode" \
         -type f \
         -name "*_${BARCODE}.mods.sorted.bam" \
         >> "${LIST}" || true

    # single libraries
    if [[ "${SAMPLE}" == "G0-C3" ]]; then
        echo "${ROOT}/G0-C3/G0-C3.mods.sorted.bam" >> "${LIST}"
    fi

    if [[ "${SAMPLE}" == "G5-C2" ]]; then
        echo "${ROOT}/G5-C2/G5-C2.mods.sorted.bam" >> "${LIST}"
    fi

    # remove missing paths just in case
    awk 'system("[ -s \"" $0 "\" ]") == 0' "${LIST}" | sort -u > "${LIST}.tmp"
    mv "${LIST}.tmp" "${LIST}"

    N=$(wc -l < "${LIST}")
    if [[ "${N}" -eq 0 ]]; then
        echo "[ERROR] No BAMs found for ${SAMPLE} ${BARCODE}" >&2
        exit 1
    fi

    echo "[merge] ${SAMPLE}: ${N} BAMs" >&2
    cat "${LIST}" >&2

    samtools merge \
        -@ "${THREADS}" \
        -f \
        -b "${LIST}" \
        "${OUTBAM}"

    samtools index -@ "${THREADS}" "${OUTBAM}"

    echo "[done] ${OUTBAM}" >&2
done
