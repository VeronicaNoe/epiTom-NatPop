#!/bin/bash
#SBATCH --output=logs/%x_%j.out         # Output file (job name and ID)
#SBATCH --error=logs/%x_%j.err          # Error file
#SBATCH --time=4-00:00:00               # Walltime (4 days)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --mem=64G			# Memory
#SBATCH --mail-type=END,FAIL		# Notifications for job done & fail
#SBATCH --mail-user=veronicanoeibanez@gmail.com
#SBATCH --job-name=preep
#SBATCH --cpus-per-task=32

# Run Bismark with paired-end reads
echo "Job started at: $(date)" >> logs/processingTime.out
make -j1
echo "Job ended at: $(date)" >> logs/processingTime.out
