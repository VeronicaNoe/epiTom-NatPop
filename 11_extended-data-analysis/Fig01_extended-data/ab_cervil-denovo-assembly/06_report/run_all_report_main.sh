#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-24}"
PROJECT="/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly"
REPORT="${PROJECT}/06_report"
PYTHON_BIN="/home/vibanez/save/miniconda3/envs/tools/bin/python"

cd "${REPORT}"

echo "[1/7] assembly summary table" >&2
"${PYTHON_BIN}" review_cervil_report_metrics.py

echo "[2/7] dotplot alignments" >&2
bash review_cervil_report_dotplots.sh

echo "[3/7] dotplot plots" >&2
"${PYTHON_BIN}" review_cervil_plot_dotplots.py

echo "[4/7] liftoff" >&2
bash review_cervil_liftoff_final.sh

echo "[5/7] liftoff summary" >&2
"${PYTHON_BIN}" review_cervil_liftoff_summary.py

echo "[6/7] initialize empty EDTA summary if missing" >&2
[[ -f edta_summary.tsv ]] || printf "source_file\tline\n" > edta_summary.tsv

echo "[7/7] markdown report" >&2
"${PYTHON_BIN}" review_cervil_make_report.py

echo "[done] main report pipeline finished" >&2
