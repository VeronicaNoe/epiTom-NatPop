#!/bin/bash
SAMPLE_DIR="/mnt/disk2/vibanez/natural-experimental_DMR/ac_output"
DMR_DIR="/mnt/disk2/vibanez/natural-experimental_DMR/aa_DMR"
CHR_SIZE="/home/vibanez/bin/chr.size.bed"
OUT_DIR="/mnt/disk2/vibanez/natural-experimental_DMR"
C_DMR="$( ls $SAMPLE_DIR/*C-DMR.bed | tr '\n' '\t' )"
CG_DMR="$( ls $SAMPLE_DIR/*CG-DMR.bed | tr '\n' '\t' )"

bedtools unionbedg -i $C_DMR -empty -g $CHR_SIZE -filler NA |\
intersectBed -a - -b $DMR_DIR/"C-DMR.bed" > $OUT_DIR/"ko-C-DMR.tmp"

bedtools unionbedg -i $CG_DMR -empty -g $CHR_SIZE -filler NA |\
intersectBed -a - -b $DMR_DIR/"CG-DMR.bed" > $OUT_DIR/"ko-CG-DMR.tmp"

ls $SAMPLE_DIR/*"C-DMR.bed" |\
sed 's./mnt/disk2/vibanez/natural-experimental_DMR/ac_output/..g' |\
sed 's/_C-DMR.bed//g'  > $OUT_DIR/00_colNames.tsv

cat $OUT_DIR/00_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }'|\
awk '{OFS="\t"}{print "chr","start","end", $0}'|\
cat - $OUT_DIR/"ko-C-DMR.tmp" | sed 's/ /\t/g' > $OUT_DIR/"ko-C-DMR.bed"

cat $OUT_DIR/00_colNames.tsv | awk 'BEGIN { ORS = "\t" } { print }'|\
awk '{OFS="\t"}{print "chr","start","end", $0}'|\
cat - $OUT_DIR/"ko-CG-DMR.tmp" | sed 's/ /\t/g' > $OUT_DIR/"ko-CG-DMR.bed"

rm $OUT_DIR/*".tmp"

