#!/bin/bash
echo '##############  ACTIVATE PLINK2 ENVIRONMENT  ################'

INDIR="/mnt/disk2/vibanez/06_get-meth-vcf/ab_epialleles"
OUTDIR="/mnt/data6/vibanez/mVCF/ag_pairwise-divergences/report"
SAMPLE_LIST="${INDIR}/toHeader-from-snps"
SAMPLES="$( cat $SAMPLE_LIST | tr '\n' '\t')"

echo "1. C-DMR differences"
echo "	1.1. get plink files"
plink2 --vcf ${INDIR}/general.C-DMR.mvcf.gz --make-bed --allow-extra-chr --out ${OUTDIR}/general.C-DMR
echo "	1.2  get differences"
plink2	--bfile ${OUTDIR}/general.C-DMR \
	--double-id \
	--allow-extra-chr\
	--threads 80 \
	--sample-diff pairwise ids=$SAMPLES \
	--out ${OUTDIR}/general.C-DMR

echo "2. CG-DMR differences"
echo "  2.1. get plink files"
plink2 --vcf ${INDIR}/general.CG-DMR.mvcf.gz --make-bed --allow-extra-chr --out ${OUTDIR}/general.CG-DMR
echo "  2.2  get differences"
plink2  --bfile ${OUTDIR}/general.CG-DMR \
        --double-id \
        --allow-extra-chr\
        --threads 80 \
        --sample-diff pairwise ids=$SAMPLES \
        --out ${OUTDIR}/general.CG-DMR

