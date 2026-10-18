#!/usr/bin/env bash
#SBATCH --job-name=cervil-report-main
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --export=NONE
#SBATCH --output=/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report/logs/report_main_%j.out
#SBATCH --error=/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report/logs/report_main_%j.err

set -euo pipefail

export HOME="/home/vibanez"
export USER="vibanez"
export LOGNAME="vibanez"
export SHELL="/bin/bash"
export PATH="/usr/local/bin:/usr/bin:/bin:${HOME}/.local/bin:${HOME}/work/bin"
export THREADS="${SLURM_CPUS_PER_TASK}"

source /etc/profile.d/modules.sh
module load bioinfo/Minimap/2-2.28

which python3
python3 --version

mkdir -p /work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report/logs
cd /work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly/06_report

bash run_all_report_main.sh
