#!/usr/bin/env bash
#SBATCH --job-name=asm-qc
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err


export THREADS="${SLURM_CPUS_PER_TASK}"

bash ~/work/bin/review_cervil-denovo-assembly_qc.sh

