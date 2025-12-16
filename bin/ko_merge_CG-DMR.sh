#!/bin/bash
CHR="$( cat $1 | cut -d'_' -f1 )"
KO="$( cat $1 | cut -d'_' -f2 )"
echo $CHR $KO
SAMPLE_DIR="09_KO-processing/09.2_biseq-processing/aa_data"
OUT_DIR="09_KO-processing/09.2_biseq-processing/ad_merge"
# intersect windows and samples
echo "------ ko-C-DMR"
cat $SAMPLE_DIR/*$CHR.$KO.bed | awk '{OFS="\t"}{print $1,$2,$3,$4}'| sortBed -i - | mergeBed -c 4,5 -o sum,sum | intersectBed -b $KO_PATH/$CHR"_"$KO.C-DMR.bed -a - -wa -wb|awk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | sortBed -i - | mergeBed -c 4,5 -o sum,sum | awk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' | sed 's/,/./g' > $OUT_DIR/$ACC"_"$CHR.$KO"-C-DMR.bed"

