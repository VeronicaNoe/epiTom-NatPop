#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-24}"
PROJECT="/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly"
REPORT="${PROJECT}/06_report"

cd "${REPORT}"

echo "[1/3] EDTA" >&2
bash review_cervil_edta_final.sh

echo "[2/3] EDTA summary" >&2
python3 review_cervil_edta_summary.py

#echo "[3/3] remake markdown report" >&2
#python3 review_cervil_make_report.py

echo "[done] EDTA pipeline finished" >&2
