#!/usr/bin/env bash
#SBATCH --job-name=dorado-mplex
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:nvidia_a100:1	#--gres=gpu:nvidia_l40s:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

module load bioinfo/Dorado/1.0.2

make -j1

