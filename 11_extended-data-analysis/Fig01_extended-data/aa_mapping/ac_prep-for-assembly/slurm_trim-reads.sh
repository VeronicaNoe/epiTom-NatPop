#!/usr/bin/env bash
#SBATCH --job-name=dorado-mplex
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

module load bioinfo/Dorado/1.0.2
module load bioinfo/samtools/1.19
module load bioinfo/Minimap/2-2.28

make -j1

