#!/bin/sh

#SBATCH -J 01_DMR-GWAS # le nom de votre job
#SBATCH -p workq  # la partition >> in this way I am asking for those nodes with higher amount of CPUs ##long for long task
#SBATCH --cpus-per-task=40
#SBATCH --mem 30000 # mémoire (en MB)
#SBATCH -o slurm-output.%N.%j.sOut # la sortie standard
#SBATCH -e slurm-error.%N.%j.sErr # la sortie erreur standard
#SBATCH --mail-type=END # notification par e-mail à la fin en cas de réussite ou erreur
#SBATCH --mail-user=veronica-noe.ibanez@universite-paris-saclay.fr # l’adresse e-mail

make -j5
