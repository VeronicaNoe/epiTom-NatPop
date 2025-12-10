#!/bin/bash
MM="$( cat $1 )"
OUT_DIR="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/ba_markers"
if [[ ${MM} == "DMR" || ${MM} == "DMR-SNP" ]]; then
	IN_DIR="/mnt/disk2/vibanez/06_get-meth-vcf/ac_vcf-metabolome"
else
	IN_DIR="/mnt/disk2/vibanez/01_raw-data/01.2_vcfiles/ab_vcf-metabolome"
fi

plink --vcf ${IN_DIR}/${MM}"_general_leaf-metabolome_LD.vcf" \
         --recode12 \
         --threads 40 \
         --allow-extra-chr \
         --output-missing-genotype 0 \
         --double-id \
         --transpose \
         --out "${OUT_DIR}/${MM}_general_leaf_LD"

if [[ ${MM} == "DMR-SNP" ]]; then
    plink --vcf ${IN_DIR}/${MM}"_general_leaf-metabolome_LD.vcf" \
             --threads 40 \
             --double-id \
             --allow-extra-chr \
             --make-bed \
             --out ${IN_DIR}/${MM}"_general_leaf-metabolome_LD"

    plink --bfile ${IN_DIR}/${MM}"_general_leaf-metabolome_LD" \
             --recode12 \
             --threads 40 \
             --double-id \
             --allow-extra-chr \
             --output-missing-genotype 0 \
             --transpose \
             --out ${OUT_DIR}/${MM}"_general_leaf_LD"
fi

~/bin/tools/EMMAX/emmax-kin-intel64 -v -d 10 ${MM}"_general_leaf_LD"
~/bin/tools/EMMAX/emmax-kin-intel64 -v -s -d 10 ${MM}"_general_leaf_LD"
