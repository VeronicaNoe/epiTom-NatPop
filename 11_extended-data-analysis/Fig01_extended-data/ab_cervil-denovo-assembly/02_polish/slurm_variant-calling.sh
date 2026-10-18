#!/usr/bin/env bash
#SBATCH --job-name=d-polish-hifiasm
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:nvidia_a100:1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
module load bioinfo/Dorado/1.0.2
module load bioinfo/samtools/1.19
export THREADS=${SLURM_CPUS_PER_TASK}
nvidia-smi

bash ~/work/bin/review_dorado_polish_hifiasm-denovo.sh
