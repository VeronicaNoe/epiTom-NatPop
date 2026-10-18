#!/usr/bin/env bash
#SBATCH --job-name=edta
#SBATCH --cpus-per-task=24
#SBATCH --mem=192G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/report_edta_%j.out
#SBATCH --error=logs/report_edta_%j.err

set -euo pipefail

export THREADS="${SLURM_CPUS_PER_TASK}"

source /etc/profile.d/modules.sh 2>/dev/null || true
module purge
module load devel/Miniconda/Miniconda3
module load bioinfo/EDTA/2.2.2


bash run_edta_only.sh
