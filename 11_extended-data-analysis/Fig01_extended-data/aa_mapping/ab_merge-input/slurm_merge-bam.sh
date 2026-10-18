#!/usr/bin/env bash
#SBATCH --job-name=merge-cervil-bam
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

module load bioinfo/samtools/1.19

bash /home/vibanez/work/bin/review_merge-bam-from-cervil.sh
