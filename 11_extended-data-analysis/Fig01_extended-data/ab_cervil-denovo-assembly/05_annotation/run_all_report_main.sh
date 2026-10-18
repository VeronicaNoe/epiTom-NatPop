#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-24}"
PROJECT="/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly"
REPORT="${PROJECT}/05_annotation"
PYTHON_BIN="/home/vibanez/save/miniconda3/envs/tools/bin/python"

cd "${REPORT}"

echo "[4/7] liftoff" >&2
bash review_cervil_liftoff_final.sh

echo "[5/7] liftoff summary" >&2
"${PYTHON_BIN}" review_cervil_liftoff_summary.py

echo "[done] main report pipeline finished" >&2
