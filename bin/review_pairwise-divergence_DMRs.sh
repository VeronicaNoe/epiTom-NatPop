#!/bin/bash
echo '##############  ACTIVATE PLINK2 ENVIRONMENT  ################'
INDIR="/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/af_DMR-epi-rate"
OUTDIR="/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/ah_pairwise-divergences"

SAMPLE_LIST="/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/aa_results/sample_names.txt"
SAMPLES="$( cat $SAMPLE_LIST | tr '\n' '\t')"

#echo "1. C-DMR differences"
echo "	1.1. get plink files"
plink2 --vcf ${INDIR}/04.1_dmr_binary_states.C-DMR.vcf.gz --make-bed --allow-extra-chr --out ${OUTDIR}/general.C-DMR-ONT

echo "	1.2  get differences"
plink2	--bfile ${OUTDIR}/general.C-DMR-ONT \
	--double-id \
	--geno 0.1 \
	--allow-extra-chr\
	--threads 80 \
	--sample-diff pairwise ids=$SAMPLES \
	--out ${OUTDIR}/general.C-DMR-ONT

echo "2. CG-DMR differences"
echo "  2.1. get plink files"
plink2 --vcf ${INDIR}/04.1_dmr_binary_states.CG-DMR.vcf.gz --make-bed --allow-extra-chr --out ${OUTDIR}/general.CG-DMR-ONT
echo "  2.2  get differences"
plink2  --bfile ${OUTDIR}/general.CG-DMR-ONT \
        --double-id \
	--geno 0.1 \
        --allow-extra-chr\
        --threads 80 \
        --sample-diff pairwise ids=$SAMPLES \
        --out ${OUTDIR}/general.CG-DMR-ONT

