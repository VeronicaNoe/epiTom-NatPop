#!/bin/bash
echo '##############  ACTIVATE PLINK2 ENVIRONMENT  ################'

VCFILE="/mnt/disk2/vibanez/01_raw-data/01.2_vcfiles/aa_download-vcf"
OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig2/ae_pairwise-divergences"
SAMPLE_LIST="${VCFILE}/sample_names.txt"
SAMPLES="$( cat $SAMPLE_LIST | tr '\n' '\t')"

#echo "SNPs differences"
#echo "  2.1. get plink files"
plink2 --vcf ${VCFILE}/SNPs_184acc-wo_indels.vcf.gz --make-bed  --allow-extra-chr --max-alleles 2 --out $OUTDIR/general.SNPs

#echo "  2.2  get differences"
plink2  --bfile ${OUTDIR}/general.SNPs \
        --double-id \
        --allow-extra-chr\
        --threads 60 \
        --sample-diff pairwise ids=$SAMPLES \
        --out $OUTDIR/general.SNPs

echo "SNPs differences"
echo "  2.1. get plink files"
#plink2 --vcf /mnt/data6/vibanez/SNPs/vcfiles/snpEff/synonymous_variants/synonymous-modifier_SNPs.vcf.gz --make-bed  --allow-extra-chr --max-alleles 2 --out $OUTDIR/synonymous-modifier_SNPs

echo "  2.2  get differences"
#plink2  --bfile $OUTDIR/synonymous-modifier_SNPs \
#        --double-id \
#        --allow-extra-chr\
#        --threads 60 \
#        --sample-diff pairwise ids=$SAMPLES \
#        --out $OUTDIR/synonymous-modifier_SNPs

