#!/usr/bin/env bash
#SBATCH --job-name=var-call
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

module load bioinfo/DeepVariant/1.6.1
module load containers/Apptainer/1.4.1
export THREADS=${SLURM_CPUS_PER_TASK}

make -j1

