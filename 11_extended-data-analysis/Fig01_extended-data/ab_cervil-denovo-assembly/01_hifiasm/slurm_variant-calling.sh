#!/usr/bin/env bash
#SBATCH --job-name=hifiasm-cervil
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load bioinfo/hifiasm
export THREADS=${SLURM_CPUS_PER_TASK}

bash ~/work/bin/review_hifiasm_denovo.sh

