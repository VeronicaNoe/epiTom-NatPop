#!/bin/bash
SAMPLE="$( cat $1 )"
INDIR="/mnt/disk2/vibanez/descriptive-stat/ac_ctxt-GWAS/sig"
OUTDIR="/mnt/disk2/vibanez/descriptive-stat/ac_ctxt-GWAS/analysis/aa_merge"
ANNODIR="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"

echo $SAMPLE
#create the bed for the DMR
cat $INDIR/$SAMPLE'.mQTL' | tr ':' '\t' | awk -v SAMPLE=$SAMPLE '{OFS="\t"}{print $1, $2, $2+1, $3, $4, $5, SAMPLE}' |\
closestBed -a - -b $ANNODIR/allGeneNames-Function.bed -t all -d |  awk '{OFS="\t"} { print $1, $2, $3, $6, $7, $11,$12,$13, $14, $15}' > $OUTDIR/$SAMPLE.annotated
