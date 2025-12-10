#!/bin/bash
SAMPLE="$( cat $1 | cut -d '_' -f1 )"
CHR="$( cat $1 | cut -d'_' -f2)"
WD="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/aa_data"
OUTPUT="$WD/tmp"
ANNO="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
CG="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/aa_data"
C="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ah_merge_C-DMR/aa_data"

echo "======== exon"
cat $CG/$SAMPLE"_leaf_Biseq_CG_"$CHR".methylation.bed" | sortBed -i - | intersectBed -a - -b $ANNO/"exon.merged.anno" -wa > $OUTPUT/$SAMPLE"_"$CHR".CG-exon"
cat $C/$SAMPLE"_leaf_Biseq_C"*$CHR".methylation.bed" | sortBed -i - | mergeBed -i - -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1, $2, $3, $4/($4+$5)}'|\
intersectBed -a - -b $ANNO/"exon.merged.anno" -wa > $OUTPUT/$SAMPLE"_"$CHR".C-exon"

cat $CG/$SAMPLE"_leaf_Biseq_CG_"$CHR".methylation.bed" | sortBed -i - | intersectBed -a - -b $ANNO/"exon-TE.merged.anno" -wa > $OUTPUT/$SAMPLE"_"$CHR".CG-exon-TE"
cat $C/$SAMPLE"_leaf_Biseq_C"*$CHR".methylation.bed" | sortBed -i - | mergeBed -i - -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1, $2, $3, $4/($4+$5)}'|\
intersectBed -a - -b $ANNO/"exon-TE.merged.anno" -wa > $OUTPUT/$SAMPLE"_"$CHR".C-exon-TE"

intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".CG-exon" -b $OUTPUT/$SAMPLE"_"$CHR".C-exon" | awk '{OFS="\t"}{print $1, $2, $3, $4, $5="teM"}' > $OUTPUT/$SAMPLE"_"$CHR"_PE-exon.bed" # if there is 1 C-DMR over the gene = teM
intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".CG-exon" -b $OUTPUT/$SAMPLE"_"$CHR"_PE-exon.bed" -v | awk '{OFS="\t"}{if ($4==0) $5="UM"; else $5="gbM"; print $1, $2, $3, $4, $5}' > $OUTPUT/$SAMPLE"_"$CHR"_only-CG-exon.bed"
intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".C-exon" -b $OUTPUT/$SAMPLE"_"$CHR"_PE-exon.bed" -v |  awk '{OFS="\t"}{if ($4==0) $5="UM"; else $5="teM"; print $1, $2, $3, $4, $5}' > $OUTPUT/$SAMPLE"_"$CHR"_only-C-exon.bed"

intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".CG-exon-TE" -b $OUTPUT/$SAMPLE"_"$CHR".C-exon-TE" | awk '{OFS="\t"}{print $1, $2, $3, $4, $5="teM"}' > $OUTPUT/$SAMPLE"_"$CHR"_PE-exon-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".CG-exon-TE" -b $OUTPUT/$SAMPLE"_"$CHR"_PE-exon-TE.bed" -v | awk '{OFS="\t"}{if ($4==0) $5="UM"; else $5="gbM"; print $1, $2, $3, $4, $5}' > $OUTPUT/$SAMPLE"_"$CHR"_only-CG-exon-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".C-exon-TE" -b $OUTPUT/$SAMPLE"_"$CHR"_PE-exon-TE.bed" -v |  awk '{OFS="\t"}{if ($4==0) $5="UM"; else $5="teM"; print $1, $2, $3, $4, $5}' > $OUTPUT/$SAMPLE"_"$CHR"_only-C-exon-TE.bed"

echo "###### intron"
cat $CG/$SAMPLE"_leaf_Biseq_CG_"$CHR".methylation.bed" | sortBed -i - | intersectBed -a - -b $ANNO/"intron.merged.anno" -wa > $OUTPUT/$SAMPLE"_"$CHR".CG-intron"
cat $C/$SAMPLE"_leaf_Biseq_C"*$CHR".methylation.bed" | sortBed -i - |mergeBed -i - -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1, $2, $3, $4/($4+$5)}'|\
intersectBed -a - -b $ANNO/"intron.merged.anno" -wa > $OUTPUT/$SAMPLE"_"$CHR".C-intron"

cat $CG/$SAMPLE"_leaf_Biseq_CG_"$CHR".methylation.bed" | sortBed -i - | intersectBed -a - -b $ANNO/"intron-TE.merged.anno" -wa > $OUTPUT/$SAMPLE"_"$CHR".CG-intron-TE"
cat $C/$SAMPLE"_leaf_Biseq_C"*$CHR".methylation.bed" | sortBed -i - |mergeBed -i - -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1, $2, $3, $4/($4+$5)}'|\
intersectBed -a - -b $ANNO/"intron-TE.merged.anno" -wa > $OUTPUT/$SAMPLE"_"$CHR".C-intron-TE"

intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".CG-intron" -b $OUTPUT/$SAMPLE"_"$CHR".C-intron" | awk '{OFS="\t"}{print $1, $2, $3, $4, $5="teM"}' > $OUTPUT/$SAMPLE"_"$CHR"_PE-intron.bed"
intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".CG-intron" -b $OUTPUT/$SAMPLE"_"$CHR"_PE-intron.bed" -v | awk '{OFS="\t"}{if ($4==0) $5="UM"; else $5="gbM"; print $1, $2, $3, $4, $5}' > $OUTPUT/$SAMPLE"_"$CHR"_only-CG-intron.bed"
intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".C-intron" -b $OUTPUT/$SAMPLE"_"$CHR"_PE-intron.bed" -v |  awk '{OFS="\t"}{if ($4==0) $5="UM"; else $5="teM"; print $1, $2, $3, $4, $5}' > $OUTPUT/$SAMPLE"_"$CHR"_only-C-intron.bed"

intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".CG-intron-TE" -b $OUTPUT/$SAMPLE"_"$CHR".C-intron-TE" | awk '{OFS="\t"}{print $1, $2, $3, $4, $5="teM"}' > $OUTPUT/$SAMPLE"_"$CHR"_PE-intron-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".CG-intron-TE" -b $OUTPUT/$SAMPLE"_"$CHR"_PE-intron-TE.bed" -v | awk '{OFS="\t"}{if ($4==0) $5="UM"; else $5="gbM"; print $1, $2, $3, $4, $5}' > $OUTPUT/$SAMPLE"_"$CHR"_only-CG-intron-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE"_"$CHR".C-intron-TE" -b $OUTPUT/$SAMPLE"_"$CHR"_PE-intron-TE.bed" -v |  awk '{OFS="\t"}{if ($4==0) $5="UM"; else $5="teM"; print $1, $2, $3, $4, $5}' > $OUTPUT/$SAMPLE"_"$CHR"_only-C-intron-TE.bed"

# epiallele state
cat $OUTPUT/$SAMPLE"_"$CHR*"-intron.bed" | sortBed -i - | awk '{OFS="\t"}{print}' > $WD/$SAMPLE"_"$CHR"_epialleles-over-intron.bed"
cat $OUTPUT/$SAMPLE"_"$CHR*"-exon.bed" | sortBed -i -  | awk '{OFS="\t"}{print}' > $WD/$SAMPLE"_"$CHR"_epialleles-over-exon.bed"
cat $OUTPUT/$SAMPLE"_"$CHR*"-intron-TE.bed" | sortBed -i - | awk '{OFS="\t"}{print}' > $WD/$SAMPLE"_"$CHR"_epialleles-over-intron-TE.bed"
cat $OUTPUT/$SAMPLE"_"$CHR*"-exon-TE.bed" | sortBed -i - | awk '{OFS="\t"}{print}' > $WD/$SAMPLE"_"$CHR"_epialleles-over-exon-TE.bed"
