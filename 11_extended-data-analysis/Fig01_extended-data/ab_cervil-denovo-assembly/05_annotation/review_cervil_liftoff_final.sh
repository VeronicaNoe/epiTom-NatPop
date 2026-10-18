#!/usr/bin/env bash
set -euo pipefail

PROJECT="/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly"
THREADS="${THREADS:-24}"

REF_FA="/home/vibanez/work/01_mapping/01.0_ref/SL5.0.fasta"
REF_GFF="/home/vibanez/work/01_mapping/01.0_ref/SL5.0.gff3"

FINAL_ASM="${PROJECT}/02_polish/cervil.polished.fasta"
LIFTOFF_BIN="/home/vibanez/.local/bin/liftoff"

OUTDIR="${PROJECT}/05_annotation/liftoff_final"
LOGDIR="${PROJECT}/05_annotation/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"

"${LIFTOFF_BIN}" \
  "${FINAL_ASM}" \
  "${REF_FA}" \
  -g "${REF_GFF}" \
  -o "${OUTDIR}/cervil.final.liftoff.gff3" \
  -u "${OUTDIR}/cervil.final.liftoff.unmapped.txt" \
  -p "${THREADS}" \
  > "${LOGDIR}/liftoff_final.stdout.log" \
  2> "${LOGDIR}/liftoff_final.stderr.log"
