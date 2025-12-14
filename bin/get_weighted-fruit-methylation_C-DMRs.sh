#!/bin/bash
SAMPLE="$( cat $1 | cut -d'_' -f1,2)"
CHR="$( cat $1 | cut -d'_' -f2 )"
SAMPLE_DIR="08_fruit-processing/08.1_DMR-classification/08.2_merge-DMRs/ab_weight-methylation/ba_data"
OUT_DIR="08_fruit-processing/08.1_DMR-classification/08.2_merge-DMRs/ab_weight-methylation"

#merge the contexts in one, sort, merge overlapping positions, get weigted methylation
cat ${SAMPLE_DIR}/${SAMPLE}*${CHR}.filtered.bed | sort -k 1.4n -k 2n | gawk '{OFS="\t"}{print $1,$2,$3,$4,$5}' |\
mergeBed -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
sed 's/,/./g'> ${OUT_DIR}/${CHR}"_"${SAMPLE}_C-DMR.weighted.methylation.bed


