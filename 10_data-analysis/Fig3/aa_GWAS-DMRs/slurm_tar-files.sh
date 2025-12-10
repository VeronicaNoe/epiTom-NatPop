#!/usr/bin/env bash
#SBATCH --job-name=compress_bb_bd
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --output=compress_bb_bd.%j.out
#SBATCH --error=compress_bb_bd.%j.err


# Option 1: two separate tarballs
tar -cvf - bb_phenotype | pigz -p 16 > bb_phenotype.tar.gz
tar -cvf - bd_results   | pigz -p 16 > bd_results.tar.gz

