#!/bin/bash
SAMPLE="$( cat $1 )"
INDIR="/mnt/disk2/vibanez/03_biseq-processing/03.4_filtering/ac_filter"
OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig1/ab_global-methylation"

cat ${INDIR}/${SAMPLE}"_"*.filtered.bed  | awk '{OFS="\t"}{print $1, $2, $2+1,$4,$5}'|\
sortBed -i - |\
mergeBed -i - -d 100000000 -c 4,5,2 -o sum,sum,count > ${OUTDIR}/$SAMPLE.weigthed-methylation-genome-wide
