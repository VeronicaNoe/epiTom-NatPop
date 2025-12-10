#!/bin/bash
SAMPLE="$( cat $1 )"
echo $SAMPLE
INDIR='/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ad_annotated-DMRs'
OTUDIR='/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ae_mergedAnnotation'
TMP='/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ae_mergedAnnotation/tmp'

ANNO="$( cat $1 | cut -d '.' -f1 )"
DMR="$( cat $1 | cut -d '.' -f2 )"

#ls $INDIR/*_$SAMPLE.methylation

cat $INDIR/*_$SAMPLE.methylation |\
awk -v ANNO="$ANNO" -v DMR="$DMR" '{OFS="\t"}{ print $1, $2, $3, $4=DMR, $5=ANNO}' >  $TMP/$SAMPLE.tmp
