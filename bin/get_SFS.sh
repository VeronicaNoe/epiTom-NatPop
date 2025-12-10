#!/bin/bash
SAMPLE="$( cat $1 )"
WD="/mnt/data6/vibanez/SNPs/vcfiles"
#INPUT="$WD/01_data" # only general SNPs_184acc-wo_indels.vcf.gz
#NAME="allAcc_general_SNPs"
INPUT="/mnt/data6/vibanez/SNPs/vcfiles/ad_annotation/bb_output"
ANN="$( cat $1 | cut -d'_' -f3 )"
NAME="allAcc_${ANN}_SNPs"
OUTPUT="$WD/ag_allele-spectrum/bb_output"
echo $SAMPLE
vcftools --gzvcf $INPUT/$SAMPLE.vcf.gz --counts --out $OUTPUT/$NAME
cat $OUTPUT/$NAME.frq.count | awk '{if(NR>1) {print}}'| awk '{if ($3<=2){print}}'|sed 's/:/\t/g' |\
awk -F '\t' 'BEGIN {OFS=FS} {if ($6 <= $8) {print $1, $2, $4, $6} else {print $1,$2,$4, $8}}' > $OUTPUT/$NAME.edited-count
#plink2 --bfile $INPUT/$SAMPLE  \
#	--allow-extra-chr \
#	--freq refbins-file=/users/bioinfo/vibanez/bin/ref-allele-bins alt1bins-file=/users/bioinfo/vibanez/bin/allele-bins \
#	--out $OUTPUT/$SAMPLE
