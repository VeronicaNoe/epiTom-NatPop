#!/usr/bin/env bash
set -euo pipefail

PROJECT="/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly"
REF="/home/vibanez/work/01_mapping/01.0_ref/SL5.0.fasta"
OUTDIR="${PROJECT}/06_report"
THREADS="${THREADS:-24}"

mkdir -p "${OUTDIR}"

declare -A ASM
ASM[hifiasm_raw]="${PROJECT}/01_hifiasm/cervil.primary.fasta"
ASM[dorado_polish]="${PROJECT}/02_polish/cervil.polished.fasta"
ASM[ragtag]="${PROJECT}/03_ragtag/ragtag.scaffold.fasta"

for STATE in hifiasm_raw dorado_polish ragtag; do
    QUERY="${ASM[$STATE]}"
    [[ -s "${QUERY}" ]] || { echo "[ERROR] missing ${QUERY}" >&2; exit 1; }

    minimap2 -x asm5 -t "${THREADS}" "${REF}" "${QUERY}" \
        > "${OUTDIR}/${STATE}.vs_SL5.paf" \
        2> "${OUTDIR}/${STATE}.vs_SL5.minimap2.log"
done
