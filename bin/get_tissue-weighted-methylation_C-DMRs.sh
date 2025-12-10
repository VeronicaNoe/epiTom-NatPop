#!/bin/bash
SAMPLE="$( cat $1 | cut -d'_' -f1,2)"
CHR="$( cat $1 | cut -d'_' -f3 )"
echo "------" $SAMPLE " " $CHR
SAMPLEPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/ae_merge_C-DMR/aa_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/ae_merge_C-DMR/bd_output"

#merge the contexts in one, sort, merge overlapping positions, get weigted methylation
cat $SAMPLEPATH/$SAMPLE*$CHR.methylation.bed | sort -k 1.4n -k 2n |\
awk '{OFS="\t"}{print $1,$2,$3,$4,$5}' | mergeBed -c 4,5 -o sum,sum |\
awk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
sed 's/,/./g' > $OUTPATH/$CHR"_"$SAMPLE.weighted.methylation.bed

