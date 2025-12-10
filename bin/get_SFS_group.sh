#!/bin/bash
SAMPLE="$( cat $1 )"
WD="/mnt/data6/vibanez/SNPs/vcfiles"
INPUT="$WD/aa_split-vcf-group/bb_output" # groups
OUTPUT="$WD/ag_allele-spectrum/bb_output"
echo $SAMPLE
vcftools --gzvcf $INPUT/$SAMPLE"-allAcc_SNPs.recode.vcf.gz" --counts --out $OUTPUT/$SAMPLE.SNPs
cat $OUTPUT/$SAMPLE.SNPs.frq.count | awk '{if(NR>1) {print}}'| awk '{if ($3<=2){print}}'|sed 's/:/\t/g' |\
awk -F '\t' 'BEGIN {OFS=FS} {if ($6 <= $8) {print $1, $2, $4, $6} else {print $1,$2,$4, $8}}' > $OUTPUT/$SAMPLE.SNPs.edited-count
#plink2 --bfile $INPUT/$SAMPLE  \
#	--allow-extra-chr \
#	--freq refbins-file=/users/bioinfo/vibanez/bin/ref-allele-bins alt1bins-file=/users/bioinfo/vibanez/bin/allele-bins \
#	--out $OUTPUT/$SAMPLE
