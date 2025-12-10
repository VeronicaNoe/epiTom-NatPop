#!/bin/bash
SAMPLE="$( cat $1 | cut -d '_' -f1)"

ALL_FEAT_DIR="/mnt/data6/vibanez/mVCF/ac_calculate-ePI/bb_output"
OUTDIR="/mnt/new_disk2/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bd_mergePI/bb_output"
TMPDIR="/mnt/new_disk2/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bd_mergePI/bb_output/tmp"

echo "1- Starting"
for file in $(ls $ALL_FEAT_DIR/*$SAMPLE"_SNPs.windowed.pi" | sed 's./mnt/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bb_general..g' ); do
cat $ALL_FEAT_DIR/$file | awk '{OFS="\t"} NR > 1 {print $1,$2,$3,$5}' > $TMPDIR/$file.tmp
done
# get all the positions
echo "2- Get positions"
cat $TMPDIR/*$SAMPLE"_SNPs.windowed.pi.tmp" | sortBed -i - | mergeBed -c 4 -o count | awk '{OFS="\t"}{print $1,$2,$3,$4}' > $TMPDIR/$SAMPLE"_allFeatures_pi-collapsed.bed"
# merge samples
echo "   "$SAMPLE
for file in $(ls $TMPDIR/*$SAMPLE"_SNPs.windowed.pi.tmp" | sed 's./mnt/new_disk2/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bd_mergePI/bb_output/tmp/..g' ); do
echo "     "$file
cat $TMPDIR/$file | intersectBed -wb -a - -b $TMPDIR/$SAMPLE"_allFeatures_pi-collapsed.bed" -nonamecheck |\
awk '{OFS="\t"}{print $1,$2,$3,$4}' > $TMPDIR/${file%%.tmp}.bed
done

echo '3- Merge samples'
CHR_SIZE="/users/bioinfo/vibanez/bin/chr.size.bed"
SAMPLE_LIST="$( ls $TMPDIR/*$SAMPLE"_SNPs.windowed.pi.bed" | tr '\n' '\t' )"
# cat $CHR_SIZE
#echo $SAMPLE_LIST
bedtools unionbedg -i $SAMPLE_LIST -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b $TMPDIR/$SAMPLE"_allFeatures_pi-collapsed.bed" > $TMPDIR/$SAMPLE"_allFeatures_merged-pi.tmp"
echo '  3.1- Adding colnames'
ls $ALL_FEAT_DIR/*$SAMPLE"_SNPs.windowed.pi" | sed 's./mnt/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bb_general/..g' |\
cut -d'_' -f1 > $OUTDIR/00_colNames_$SAMPLE.tsv

cat $OUTDIR/00_colNames_$SAMPLE.tsv | awk 'BEGIN { ORS = "\t" } { print }'| awk '{OFS="\t"}{print "chr","start","end", $0}'|\
cat - $TMPDIR/$SAMPLE"_allFeatures_merged-pi.tmp"  > $OUTDIR/"SNPs_general_"$SAMPLE"_merged-pi.bed"

rm $TMPDIR/$SAMPLE"_allFeatures_merged-pi.tmp"
