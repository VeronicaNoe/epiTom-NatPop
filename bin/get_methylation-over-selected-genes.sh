#!/bin/bash
SAMPLE_NAME="$( cat $1 | cut -d'_' -f1)"
TYPE="$( cat $1 | cut -d'_' -f2)"

OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig4/ac_get-methylation-levels-over-selected-genes"
GENELIST="/mnt/disk2/vibanez/10_data-analysis/Fig4/results/geneList_coordinates.bed"
if [[ "$TYPE" == "KO" ]]; then
  INDIR="/mnt/disk2/vibanez/09_KO-processing/09.0_biseq/ad_filtering/bb_ctxt-splitting"
else
  INDIR="/mnt/disk2/vibanez/03_biseq-processing/03.4_filtering/ac_filter"
fi
cat $INDIR/${SAMPLE_NAME}_leaf*bed | awk '{OFS="\t"}{print $1, $2, $2+1, $4, $5}' |\
intersectBed -a - -b $GENELIST -wo | sortBed -i - |\
mergeBed -i - -c 4,5,10 -o sum,sum,distinct > $OUTDIR/${SAMPLE_NAME}_methylation-over-selected-genes.bed
