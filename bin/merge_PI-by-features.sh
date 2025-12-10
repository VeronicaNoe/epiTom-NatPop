#!/bin/bash
SAMPLE="$( cat $1 | cut -d '_' -f1)"
ANNO="$( cat $1 | cut -d '_' -f2 )"

BY_FEAT_DIR="/mnt/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bc_annotated/bb_output"
OUTDIR="/mnt/new_disk2/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bd_mergePI/bc_output"
TMPDIR="/mnt/new_disk2/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bd_mergePI/bc_output/tmp"

echo "1- Starting"
for file in $(ls $BY_FEAT_DIR/"SNPs_"*"_"$ANNO"_"$SAMPLE".windowed.pi" | sed 's./mnt/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bc_annotated/bb_output/..g' ); do
ls -lh $BY_FEAT_DIR/$file
cat $BY_FEAT_DIR/$file | awk '{OFS="\t"} NR > 1 {print $1,$2,$3,$5}' > $TMPDIR/$file.tmp
done

# get all the positions
echo "2- Get positions"
cat $TMPDIR/"SNPs_"*"_"$ANNO"_"$SAMPLE".windowed.pi.tmp" | sortBed -i - | mergeBed -c 4 -o count | awk '{OFS="\t"}{print $1,$2,$3,$4}' > $TMPDIR/"SNPs_"$ANNO"_"$SAMPLE"_pi-collapsed.bed"
# merge samples
echo "   "$SAMPLE

for file in $(ls $TMPDIR/"SNPs_"*"_"$ANNO"_"$SAMPLE".windowed.pi.tmp" | sed 's./mnt/new_disk2/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bd_mergePI/bc_output/tmp/..g' ); do
echo "     "$file
cat $TMPDIR/$file | intersectBed -wb -a - -b $TMPDIR/"SNPs_"$ANNO"_"$SAMPLE"_pi-collapsed.bed" -nonamecheck |\
awk '{OFS="\t"}{print $1,$2,$3,$4}' > $TMPDIR/${file%%.tmp}.bed
done

echo '3- Merge samples'
CHR_SIZE="/users/bioinfo/vibanez/bin/chr.size.bed"
SAMPLE_LIST="$( ls $TMPDIR/"SNPs_"*"_"$ANNO"_"$SAMPLE".windowed.pi.bed" | tr '\n' '\t' )"
#cat $CHR_SIZE
echo $SAMPLE_LIST
bedtools unionbedg -i $SAMPLE_LIST -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b $TMPDIR/"SNPs_"$ANNO"_"$SAMPLE"_pi-collapsed.bed" > $TMPDIR/"SNPs_"$ANNO"_"$SAMPLE"_merged-pi.tmp"

echo '  3.1- Adding colnames'
ls $BY_FEAT_DIR/"SNPs_"*"_"$ANNO"_"$SAMPLE".windowed.pi" | sed 's./mnt/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bc_annotated/bb_output/..g' |\
cut -d'_' -f2 > $OUTDIR/00_colNames_$SAMPLE"_"$ANNO.tsv

cat $OUTDIR/00_colNames_$SAMPLE"_"$ANNO.tsv | awk 'BEGIN { ORS = "\t" } { print }'| awk '{OFS="\t"}{print "chr","start","end", $0}'|\
cat - $TMPDIR/"SNPs_"$ANNO"_"$SAMPLE"_merged-pi.tmp"  > $OUTDIR/"SNPs_"$ANNO"_"$SAMPLE"_merged-pi.bed"

rm $TMPDIR/"SNPs_"$ANNO"_"$SAMPLE"_merged-pi.tmp"
