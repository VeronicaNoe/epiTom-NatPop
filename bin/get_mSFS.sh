#!/bin/bash
SAMPLE="$( cat $1 )"
#INPUT="/mnt/data6/vibanez/mVCF/ab_split-mvcf-group/aa_data" # for data split by groups
INPUT="/mnt/data6/vibanez/mVCF/aa_data" # for general
OUTPUT="/mnt/data6/vibanez/SNPs/vcfiles/ag_allele-spectrum/bb_output"
NAME="allAcc_${SAMPLE}"
ls $INPUT/$SAMPLE
vcftools --gzvcf $INPUT/$SAMPLE".vcf.gz" --counts --out $OUTPUT/$NAME
cat $OUTPUT/$NAME".frq.count" | awk '{if(NR>1) {print}}'|\
awk '{if ($3<=2){print}}'|sed 's/:/\t/g' |\
awk -F '\t' 'BEGIN {OFS=FS} {if ($6 <= $8) {print $1, $2, $4, $6} else {print $1,$2,$4, $8}}' > $OUTPUT/$NAME".edited-count"
