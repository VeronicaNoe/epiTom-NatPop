#!/bin/bash
ACC="$( cat $1 | cut -d'_' -f1 )"
CHR="$( cat $1 | cut -d'_' -f2 )"
KO="$( cat $1 | cut -d'_' -f3 )"
echo $ACC $CHR $KO

SAMPLE_DIR="09_KO-processing/09.0_biseq/ad_filtering/bc_filter-window"
DMR_DIR="09_KO-processing/09.0_biseq/ae_DMR"
OUT_DIR="09_KO-processing/09.2_biseq-processing/aa_data"

# intersect windows and samples
echo "------ ko-C-DMR"
cat $SAMPLE_DIR/$ACC"_CH"*$CHR.filtered.bed | awk '{OFS="\t"}{print $1,$2,$2,$4,$5, $6}'|\
sortBed -i - | mergeBed -c 4,5 -o sum,sum |\
intersectBed -b $DMR_DIR/$CHR"_"$KO.C-DMR.bed -a - -wa -wb|awk '{OFS="\t"}{print $1,$7,$8,$4,$5}' |\
sortBed -i - | mergeBed -c 4,5 -o sum,sum | awk '{OFS="\t"}{print $1,$2,$3,$4,$5}' > $OUT_DIR/$ACC"_"$CHR.$KO"-C-DMR.bed"
echo "------ ko-CG-DMR"
cat $SAMPLE_DIR/$ACC"_CG"*$CHR.filtered.bed | awk '{OFS="\t"}{print $1,$2,$2,$4,$5, $6}'|\
bedtools subtract -a - -b $OUT_DIR/$ACC"_"$CHR.$KO"-C-DMR.bed" | intersectBed -b $DMR_DIR/$CHR"_"$KO.CG-DMR.bed -a - -wa -wb |\
awk '{OFS="\t"}{print $1,$8,$9,$4,$5}' | mergeBed -c 4,5 -o sum,sum |\
awk '{OFS="\t"}{print $1,$2,$3,$4,$5}'  > $OUT_DIR/$ACC"_"$CHR.$KO"-CG-DMR.bed"

