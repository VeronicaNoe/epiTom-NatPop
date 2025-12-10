#!/bin/bash
CHR="$( cat $1 | cut -d'_' -f1 )"
KO="$( cat $1 | cut -d'_' -f2 )"
echo $CHR $KO

#bash
#conda activate bedtools
# --- Define samples to process
SAMPLE_PATH="/mnt/disk2/vibanez/02_methylkit/ad_KO-DMR/bb_output"
#KO_PATH="/mnt/disk2/vibanez/02_methylkit/ad_KO-DMR/aa_data"
OUT="/mnt/disk2/vibanez/02_methylkit/ad_KO-DMR/bd_merge"
ls $SAMPLE_PATH/*$CHR.$KO.bed

# intersect windows and samples
#echo "------ ko-C-DMR"
cat $SAMPLE_PATH/*$CHR.$KO.bed | awk '{OFS="\t"}{print $1,$2,$3,$4}'| sortBed -i - | mergeBed -c 4,5 -o sum,sum | intersectBed -b $KO_PATH/$CHR"_"$KO.C-DMR.bed -a - -wa -wb|awk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | sortBed -i - | mergeBed -c 4,5 -o sum,sum | awk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' | sed 's/,/./g' > $OUT/$ACC"_"$CHR.$KO"-C-DMR.bed"
#echo "------ ko-CG-DMR"
#cat $SAMPLE_PATH/$ACC"_CG"*$CHR.filtered.bed | awk '{OFS="\t"}{print $1,$2,$2,$4,$5, $6}'| bedtools subtract -a - -b $OUT/$ACC"_"$CHR.$KO"-C-DMR.bed" | intersectBed -b $KO_PATH/$CHR"_"$KO.CG-DMR.bed -a - -wa -wb | awk '{OFS="\t"}{print $1,$8,$9,$4,$5}' | mergeBed -c 4,5 -o sum,sum | awk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' | sed 's/,/./g'  > $OUT/$ACC"_"$CHR.$KO"-CG-DMR.bed"

