#!/bin/bash
SAMPLES_FILE=$1
#conda activate bedtools
# --- Define samples to process
SAMPLE="$( cat "$SAMPLES_FILE" )"
SAMPLE_PATH="/mnt/disk2/vibanez/01_bismark_analysis/04_processingFiles/ac_filter/bb_output"
WINDOW="/mnt/disk2/vibanez/02_methylkit/ai_meth-description"
echo "------ bedtools format"

cat $SAMPLE_PATH/$SAMPLE"_ch"*.filtered.bed  | awk '{OFS="\t"}{print $1, $2, $2,$4,$5}'|\
sortBed -i - | intersectBed -a - -b $WINDOW/genome.10mb.window -wao|\
awk '{OFS="\t"}{print $1, $7, $8, $4, $5}' |\
sortBed -i - | mergeBed -i - -c 4,5 -o sum,sum |\
awk '{OFS="\t"}{print $1, $2, $3, $4/($4+$5)}'|\
sed 's/,/./g' > $SAMPLE.meth-10mb.windows
