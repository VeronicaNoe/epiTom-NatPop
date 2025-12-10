#!/bin/bash
TISSUE="$( cat $1 )"
#SNP="/mnt/disk2/vibanez/01_raw-data/01.2_vcfiles/aa_download-vcf/SNPs_SL2.5_185-samples_biallelic_wo_indels.vcf.gz"
SNP="/mnt/disk6/vibanez/SNPs"
OUT_DIR="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome"
zcat ${SNP}_SL2.5_185-samples_biallelic_wo_indels_with_id.vcf.gz | sed 's/SL2.50ch//g' > SNPs_leaf_w_annotation_chr_edited.vcf
bgzip SNPs_leaf_w_annotation_chr_edited.vcf

#rm tmp_SNPS_${TISSUE}_*
##
plink1.9 --vcf SNPs_leaf_w_annotation_chr_edited.vcf.gz --threads 40 --double-id --allow-extra-chr --maf 0.05 --indep-pairwise 50 5 0.5 --make-bed --out SNP_general_leaf
plink1.9 --bfile SNP_general_leaf --extract SNP_general_leaf.prune.in --make-bed --allow-extra-chr --recode vcf-fid --out SNP_general_leaf_LD

