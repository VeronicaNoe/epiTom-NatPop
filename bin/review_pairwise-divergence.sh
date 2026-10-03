#!/bin/bash
echo '##############  ACTIVATE PLINK2 ENVIRONMENT  ################'

VCFILE="/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/aa_results"
OUTDIR="/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/ah_pairwise-divergences"
SAMPLE_LIST="${VCFILE}/sample_names.txt"
SAMPLES="$( cat $SAMPLE_LIST | tr '\n' '\t')"

#echo "SNPs differences"
#echo "  2.1. get plink files"
plink2 --vcf ${VCFILE}/cervil-SNPs-ONT.GQ20_DP5.tags.vcf.gz --make-bed  --allow-extra-chr --max-alleles 2 --out $OUTDIR/general.SNPs-ONT

#echo "  2.2  get differences"
plink2  --bfile ${OUTDIR}/general.SNPs-ONT \
        --double-id \
	--geno 0.1 \
        --allow-extra-chr\
        --threads 60 \
        --sample-diff pairwise ids=$SAMPLES \
        --out $OUTDIR/general.SNPs-ONT


