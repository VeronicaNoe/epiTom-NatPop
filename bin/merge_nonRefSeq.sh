#!/bin/bash
INDIR="/mnt/disk2/vibanez/02_methylkit/nonRefSeq/aa_data"
OUTDIR="/mnt/disk2/vibanez/02_methylkit/nonRefSeq/ab_merged"
CTXT="$( echo $1 )"
CG="$( ls $INDIR/*"_leaf_Biseq_"$CTXT"_chrNonRef" | tr '\n' '\t' )"

echo "Union all bed files"
#echo $CG
bedtools unionbedg -i $CG -empty -g $INDIR/"nonRefSeq.size" -filler NA  > $OUTDIR/$CTXT"_chrNonRef.tmp"

#echo "Geting the samples Names"
#ls $INDIR/*"CG_chrNonRef" |sed 's./mnt/disk2/vibanez/02_methylkit/nonRefSeq/aa_data/..g' | sed 's/_CG_chrNonRef//g' |\
#cut -d'_' -f1,2 > $OUTDIR/00_colNames.tsv

echo "Adding samples Names"
cat $OUTDIR/00_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }'| awk '{OFS="\t"}{print "chr","start","end", $0}'|\
cat - $OUTDIR/$CTXT"_chrNonRef.tmp" > $OUTDIR/$CTXT"_chrNonRef.bed"
