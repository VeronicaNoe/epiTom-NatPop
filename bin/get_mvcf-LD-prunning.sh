#!/bin/bash
DIR="/mnt/disk2/vibanez/06_get-meth-vcf/ab_epialleles"
#use plink1.9
for DMR in CG-DMR C-DMR; do
  echo " Processing $DMR..."

  # Step 1: Create initial PLINK binary files and perform LD pruning
  plink --vcf ${DIR}/${DMR}_general_leaf-metabolome.vcf \
	--threads 40 \
	--const-fid \
	--allow-extra-chr \
	--maf 0.05 \
	--make-bed \
	--indep-pairwise 100 100 0.2 \
	--out ${DIR}/${DMR}_general_leaf-metabolome

  # Step 2: Generate LD-pruned VCF
  plink --bfile ${DIR}/${DMR}_general_leaf-metabolome \
           --extract ${DIR}/${DMR}_general_leaf-metabolome.prune.in \
           --make-bed \
           --allow-extra-chr \
           --recode vcf \
           --const-fid \
           --out ${DIR}/${DMR}_general_leaf-metabolome_LD
done
