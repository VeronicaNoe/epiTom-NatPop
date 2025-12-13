#!/bin/bash
TISSUE="$( cat $1 )"
SNP="/01_raw-data/01.2_vcfiles/aa_download-vcf/SNPs_SL2.5_185-samples_biallelic_wo_indels_with_id.vcf.gz"
OUT_DIR="10_data-analysis/Fig5/aa_GWAS-metabolome"
bcftools view -S samples_${TISSUE}_metabolome -Oz -o tmp_SNPs_${TISSUE}_metabolome.vcf.gz ${SNP}
bcftools annotate --threads 40 --set-id '+%CHROM:%POS' tmp_SNPs_${TISSUE}_metabolome.vcf.gz -O z -o tmp_SNPs_${TISSUE}_metabolome_w_annotation.vcf.gz
zcat tmp_SNPs_${TISSUE}_metabolome_w_annotation.vcf.gz | sed 's/SL2.50ch//g' > SNPs_${TISSUE}_metabolome_w_annotation_chr_edited.vcf
bgzip SNPs_${TISSUE}_metabolome_w_annotation_chr_edited.vcf

rm tmp_SNPS_${TISSUE}_*
##
plink1.9 --vcf SNPs_${TISSUE}_metabolome_w_annotation_chr_edited.vcf.gz --threads 40 --double-id --allow-extra-chr --maf 0.05 --indep-pairwise 50 5 0.5 --make-bed --out SNP_general_${TISSUE}-metabolome
plink1.9 --bfile SNP_general_${TISSUE}-metabolome --extract SNP_general_${TISSUE}-metabolome.prune.in --make-bed --allow-extra-chr --recode vcf-fid --out SNP_general_${TISSUE}-metabolome_LD
plink1.9 --vcf SNP_general_${TISSUE}-metabolome_LD.vcf --recode12 --threads 40 --allow-extra-chr --output-missing-genotype 0 --const-fid --transpose --out ${OUT_DIR}/${TISSUE}/SNPs/data/SNP_general_${TISSUE}_LD

