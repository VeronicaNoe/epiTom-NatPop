#!/usr/bin/env bash
#SBATCH --job-name=sniffles
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

module load devel/Miniforge/Miniforge3
eval "$(conda shell.bash hook)"
conda activate Sniffles-v2.6.1_env

make -j1
