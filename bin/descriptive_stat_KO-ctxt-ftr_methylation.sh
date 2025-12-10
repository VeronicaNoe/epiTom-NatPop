#!/bin/bash
echo 'Activate gawk environment'
CTXT="$( cat $1 | cut -d'_' -f1 )"
KO="$( cat $1 | cut -d'_' -f2 )"
KO2="$( cat $1 | cut -d'_' -f3)"
ANNO="$( cat $1 | cut -d'.' -f1 | cut -d'_' -f4,5 )"

KODIR="/mnt/disk2/vibanez/04_natural-experimental_DMR/aa_KO"
ANNODIR="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTDIR="/mnt/disk3/vibanez/descriptive-stat/ab_annotated-bed/ko_meth"

cat $KODIR/$KO".ctxt-merged"  | gawk -v CTXT=$CTXT '{OFS="\t"}{ if ($10==CTXT) {print}}'  |\
gawk '{OFS="\t"}{ if ($13<0) {print}}' |\
intersectBed -a - -b $ANNODIR/$ANNO.anno -wb |\
intersectBed -a - -b $KODIR/$KO2".ctxt-merged" |\
gawk '{sum4 += $4; sum5 += $5; sum7 += $7; sum8 += $8}  END {print sum4, sum5, sum5/sum4, sum7, sum8, sum8/sum7}' OFS="\t" > $OUTDIR/$CTXT"_"$KO"_"$KO2"_"$ANNO".methylation"

