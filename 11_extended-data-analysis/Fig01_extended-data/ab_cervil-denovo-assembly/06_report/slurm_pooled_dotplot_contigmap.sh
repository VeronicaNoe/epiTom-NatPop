#!/usr/bin/env bash
#SBATCH --job-name=cervil-pooled-qc
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --output=/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report/logs/pooled_qc_%j.out
#SBATCH --error=/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report/logs/pooled_qc_%j.err

set -euo pipefail

export HOME="/home/vibanez"
export USER="vibanez"
export LOGNAME="vibanez"
export SHELL="/bin/bash"
export PATH="/usr/local/bin:/usr/bin:/bin:${HOME}/.local/bin:${HOME}/work/bin"

cd /work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report

mkdir -p logs

which python3
python3 --version

python3 review_cervil_pooled_dotplot_and_contigmap.py
