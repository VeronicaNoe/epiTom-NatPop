#!/bin/bash
GROUP="$( cat $1 | cut -d'_' -f1 )"
IN_DIR="01_raw-data/01.2_vcfiles/aa_download-vcf"
IN_FILE="01_raw-data/01.2_vcfiles/aa_download-vcf/SNPs_SL2.5_185-samples_biallelic_wo_indels_with_id.vcf.gz"
OUTDIR="01_raw-data/01.2_vcfiles/ac_processed-vcf"
echo $GROUP

vcftools --gzvcf $IN_FILE \
	--keep $IN_DIR/$GROUP".group" \
	--recode \
	--recode-INFO-all \
	--out $OUTDIR/$GROUP"_allAcc"

vcftools --gzvcf $IN_FILE \
        --keep $IN_DIR/$GROUP".20acc.group" \
        --recode \
        --recode-INFO-all \
        --out $OUTDIR/$GROUP"_20acc"
