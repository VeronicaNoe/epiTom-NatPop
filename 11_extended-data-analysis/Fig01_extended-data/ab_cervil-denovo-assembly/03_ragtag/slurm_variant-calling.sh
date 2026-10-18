#!/usr/bin/env bash
#SBATCH --job-name=ragtag-hifiasm
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
module load bioinfo/RagTag || true
export THREADS=${SLURM_CPUS_PER_TASK}

#make -j1
bash ~/work/bin/review_ragtag_hifiasm_denovo.sh
