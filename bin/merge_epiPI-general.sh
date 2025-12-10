#!/bin/bash
SAMPLE="$( cat $1 | cut -d '_' -f2)"
CTXT="$( cat $1 | cut -d '_' -f1)"

DIR="/mnt/data6/vibanez/mVCF/ac_calculate-ePI/ba_data"
OUTDIR="/mnt/data6/vibanez/mVCF/ac_calculate-ePI/bc_output"
TMPDIR="/mnt/data6/vibanez/mVCF/ac_calculate-ePI/tmp"

echo "1- Starting"
#remove header
for file in $(ls $DIR/*"_"$CTXT"_general_"$SAMPLE".recode.windowed.pi" | sed 's./mnt/data6/vibanez/mVCF/ac_calculate-ePI/ba_data/..g' ); do
ls -lh $DIR/$file
cat $DIR/$file | awk '{OFS="\t"} NR > 1 {print $1,$2,$3,$5}' > $TMPDIR/$file.tmp
done

# get all the positions
echo "2- Get positions"
cat $TMPDIR/*"_"$CTXT"_general_"$SAMPLE".recode.windowed.pi.tmp" | sortBed -i - | mergeBed -c 4 -o count | awk '{OFS="\t"}{print $1,$2,$3,$4}' > $TMPDIR/$CTXT"_general_"$SAMPLE"_pi-collapsed.bed"
# merge samples
echo "   "$SAMPLE
for file in $(ls $TMPDIR/*"_"$CTXT"_general_"$SAMPLE".recode.windowed.pi.tmp" | sed 's./mnt/data6/vibanez/mVCF/ac_calculate-ePI/tmp/..g' ); do
echo "     "$file
cat $TMPDIR/$file | intersectBed -wb -a - -b $TMPDIR/$CTXT"_general_"$SAMPLE"_pi-collapsed.bed" -nonamecheck |\
awk '{OFS="\t"}{print $1,$2,$3,$4}' > $TMPDIR/${file%%.tmp}.bed
done

echo '3- Merge samples'
CHR_SIZE="/users/bioinfo/vibanez/bin/chr.size.bed"
SAMPLE_LIST="$( ls $TMPDIR/*"_"$CTXT"_general_"$SAMPLE".recode.windowed.pi.bed" | tr '\n' '\t' )"
#cat $CHR_SIZE
echo $SAMPLE_LIST
bedtools unionbedg -i $SAMPLE_LIST -empty -g $CHR_SIZE -filler NA | intersectBed -a - -b $TMPDIR/$CTXT"_general_"$SAMPLE"_pi-collapsed.bed" > $TMPDIR/$CTXT"_general_"$SAMPLE"_merged-pi.tmp"

echo '  3.1- Adding colnames'
ls $DIR/*"_"$CTXT"_general_"$SAMPLE".windowed.pi" | sed 's./mnt/data6/vibanez/mVCF/ac_calculate-ePI/ba_data/..g' |\
cut -d'_' -f2 > $OUTDIR/00_colNames_$SAMPLE"_general.tsv"

cat $OUTDIR/00_colNames_$SAMPLE"_general.tsv" | awk 'BEGIN { ORS = "\t" } { print }'| awk '{OFS="\t"}{print "chr","start","end", $0}'|\
cat - $TMPDIR/$CTXT"_general_"$SAMPLE"_merged-pi.tmp"  > $OUTDIR/$CTXT"_general_"$SAMPLE"_merged-pi.bed"

rm $TMPDIR/$CTXT"_general_"$SAMPLE"_merged-pi.tmp"
