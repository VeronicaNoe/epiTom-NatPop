#!/usr/bin/env bash
#SBATCH --job-name=cervil-plotsr
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report/logs/plotsr_%j.out
#SBATCH --error=/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report/logs/plotsr_%j.err

set -euo pipefail

PROJECT="/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly"
REPORT="${PROJECT}/06_report"
OUTDIR="${REPORT}/syri_plotsr_chr12"

cd "${OUTDIR}"

source /etc/profile.d/modules.sh 2>/dev/null || true
module purge
module load devel/Miniconda/Miniconda3
module load bioinfo/plotsr/1.1.1

echo "[check tools]" >&2
which plotsr >&2
plotsr --help | head >&2 || true

if [[ ! -s cervil_vs_SL5.syri.out ]]; then
    echo "[ERROR] Missing SyRI output: ${OUTDIR}/cervil_vs_SL5.syri.out" >&2
    ls -lh "${OUTDIR}" >&2
    exit 1
fi

if [[ ! -s genomes.txt ]]; then
    echo "[ERROR] Missing genomes file: ${OUTDIR}/genomes.txt" >&2
    ls -lh "${OUTDIR}" >&2
    exit 1
fi

echo "[plotsr]" >&2
plotsr \
  --sr cervil_vs_SL5.syri.out \
  --genomes genomes.txt \
  -o "${REPORT}/cervil_SL5_syri_plotsr.pdf" \
  > plotsr.stdout.log \
  2> plotsr.stderr.log

echo "[done] plotsr finished" >&2
ls -lh "${REPORT}/cervil_SL5_syri_plotsr.pdf" >&2
