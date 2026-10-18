#!/usr/bin/env bash
set -euo pipefail

PROJECT="/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly"
THREADS="${THREADS:-24}"

FINAL_ASM="${PROJECT}/02_polish/cervil.polished.fasta"
OUTDIR="${PROJECT}/05_annotation/edta_final"
LOGDIR="${PROJECT}/05_annotation/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"
cd "${OUTDIR}"

EDTA.pl \
  --genome "${FINAL_ASM}" \
  --species others \
  --step all \
  --threads "${THREADS}" \
  > "${LOGDIR}/edta_final.stdout.log" \
  2> "${LOGDIR}/edta_final.stderr.log"
