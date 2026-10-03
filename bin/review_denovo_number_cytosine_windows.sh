#!/bin/bash
echo 'conda activate bedtools'
# --- Define samples to process
SAMPLE="$( cat "$1")"
echo $SAMPLE
SAMPLE_DIR="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/af_chr-split"
TMP="/mnt/disk1/vibanez/tmp"
mkdir -p ${TMP}
echo "------ bedtools format"
cat ${SAMPLE_DIR}/$SAMPLE.bed  | gawk '{OFS="\t"}{print $1, $2, $2+1, $3,$4,$5,$6,$7}' > $TMP/$SAMPLE.tmp
# intersect windows and samples
echo "------ bedtools by windows"
intersectBed -c -a /mnt/ssd123/vibanez/19_ont-mapping/ab_run/cervil.polished.contigs_used_in_ragtag.windows_100bp.bed -b $TMP/$SAMPLE.tmp > $TMP/$SAMPLE.counts_windows.tmp
# keep only windows with more than 3 c per strand
echo "------ geting only more than 3"
cat $TMP/$SAMPLE.counts_windows.tmp | gawk '{if ($4 >=3) print}' | sed 's/ /\t/g' > $TMP/$SAMPLE.counts_windows_filtered.tmp
# get the samples filtered coordinates to keep
echo "------ filtering"
intersectBed -a $TMP/$SAMPLE.tmp -b $TMP/$SAMPLE.counts_windows_filtered.tmp > $TMP/$SAMPLE.filtered.bed
# save in the correct format
echo "------ Saving"
cat $TMP/$SAMPLE.filtered.bed | gawk '{OFS="\t"}{print $1, $2, $4, $5, $6, $7, $8}'  > $SAMPLE.filtered.bed
# remove tmp filesa
#rm $TMP/$SAMPLE*.tmp
