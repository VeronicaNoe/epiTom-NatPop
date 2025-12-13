#!/bin/bash
# --- Define samples to process
SAMPLE="$( cat "$1")"
CHR="$( cat "$1" | grep -o 'SL2\.50ch[0-9]\{2\}' | sed 's/SL2.50//g' )"
echo $CHR
echo $SAMPLE
SAMPLE_DIR="03_biseq-processing/03.4_filtering/ab_chr-split"
TMP="tmp"
echo "------ bedtools format"
cat ${SAMPLE_DIR}/$SAMPLE.bed  | gawk '{OFS="\t"}{print $1, $2, $2+1, $3,$4,$5,$6,$7}' > $TMP/$SAMPLE.tmp
# intersect windows and samples
echo "------ bedtools by windows"
intersectBed -c -a ~/bin/tmp/$CHR.window -b $TMP/$SAMPLE.tmp > $TMP/$SAMPLE.counts_windows.tmp
# keep only windows with more than 3 c per strand
echo "------ geting only more than 3"
cat $TMP/$SAMPLE.counts_windows.tmp | gawk '{if ($4 >=3) print}' | sed 's/ /\t/g' > $TMP/$SAMPLE.counts_windows_filtered.tmp
# get the samples filtered coordinates to keep
echo "------ filtering"
intersectBed -a $TMP/$SAMPLE.tmp -b $TMP/$SAMPLE.counts_windows_filtered.tmp > $TMP/$SAMPLE.filtered.bed
# save in the correct format
echo "------ Saving"
cat $TMP/$SAMPLE.filtered.bed | gawk '{OFS="\t"}{print $1, $2, $4, $5, $6, $7, $8}'  > $SAMPLE.filtered.bed
# remove tmp files
rm $TMP/$SAMPLE*.tmp
