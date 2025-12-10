#!/bin/bash
SAMPLE="$( cat $1 | cut -d '_' -f1,2 )"
echo $CHR
WD="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/aa_data"
OUTPUT="$WD/tmp"
ANNO="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
CG="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ag_merge_CG-DMR/bb_output"
C="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ah_merge_C-DMR/bb_output"
OUTMETH="/mnt/disk2/vibanez/02_methylkit/ag_poliepiallelic-genes/methylation/aa_data"

## exon
# get exon coordinates with  complete DMRs
cat $CG/$CHR"_CG-DMR.merged.methylation.bed" | awk '{OFS="\t"} NR>1{print}' | intersectBed -b -  -a $ANNO/"exon-wo-TE.merged.anno" -wao |\
awk '{if ($7>=99) {print}}'| mergeBed -i -  > $OUTPUT/$CHR".CG-exon"
cat $C/$CHR"_C-DMR.merged.methylation.bed" | awk '{OFS="\t"} NR>1{print}' | intersectBed -b -  -a $ANNO/"exon-wo-TE.merged.anno" -wao |\
awk '{if ($7>=99) {print}}'| mergeBed -i -  > $OUTPUT/$CHR".C-exon"
cat $CG/$CHR"_CG-DMR.merged.methylation.bed" | awk '{OFS="\t"} NR>1{print}' | intersectBed -b -  -a $ANNO/"exon-TE.merged.anno" -wao |\
awk '{if ($7>=99) {print}}'| mergeBed -i -  > $OUTPUT/$CHR".CG-exon-TE"
cat $C/$CHR"_C-DMR.merged.methylation.bed" | awk '{OFS="\t"} NR>1{print}' | intersectBed -b -  -a $ANNO/"exon-TE.merged.anno" -wao |\
awk '{if ($7>=99) {print}}'| mergeBed -i -  > $OUTPUT/$CHR".C-exon-TE"

intersectBed -a $OUTPUT/$CHR".CG-exon" -b $OUTPUT/$CHR".C-exon" | awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}' > $OUTPUT/$CHR"_PE-exon.bed"
intersectBed -a $OUTPUT/$CHR".CG-exon" -b $OUTPUT/$CHR"_PE-exon.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4="gbM"}'> $OUTPUT/$CHR"_only-CG-exon.bed"
intersectBed -a $OUTPUT/$CHR".C-exon" -b $OUTPUT/$CHR"_PE-exon.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}'> $OUTPUT/$CHR"_only-C-exon.bed"
intersectBed -a $OUTPUT/$CHR".CG-exon-TE" -b $OUTPUT/$CHR".C-exon-TE" | awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}' > $OUTPUT/$CHR"_PE-exon-TE.bed"
intersectBed -a $OUTPUT/$CHR".CG-exon-TE" -b $OUTPUT/$CHR"_PE-exon-TE.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4="gbM"}'> $OUTPUT/$SAMPLE"_only-CG-exon-TE.bed"
intersectBed -a $OUTPUT/$CHR".C-exon-TE" -b $OUTPUT/$CHR"_PE-exon-TE.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}'> $OUTPUT/$SAMPLE"_only-C-exon-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-exon" -b $OUTPUT/$SAMPLE".C-exon" | awk '{OFS="\t"}{print $1, $2, $3, $4}' > $OUTMETH/$SAMPLE"_PE-exon.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-exon" -b $OUTPUT/$SAMPLE"_PE-exon.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/$SAMPLE"_only-CG-exon.bed"
intersectBed -a $OUTPUT/$SAMPLE".C-exon" -b $OUTPUT/$SAMPLE"_PE-exon.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/$SAMPLE"_only-C-exon.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-exon-TE" -b $OUTPUT/$SAMPLE".C-exon-TE" | awk '{OFS="\t"}{print $1, $2, $3, $4}' > $OUTMETH/$SAMPLE"_PE-exon-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-exon-TE" -b $OUTPUT/$SAMPLE"_PE-exon-TE.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/$SAMPLE"_only-CG-exon-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE".C-exon-TE" -b $OUTPUT/$SAMPLE"_PE-exon-TE.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/$SAMPLE"_only-C-exon-TE.bed"

## intron
intersectBed -b $CG/$SAMPLE"_leaf.CG-DMR.bed"  -a $ANNO/"intron-wo-TE.merged.anno" -wao | awk '{if ($8>=99){print}}'| mergeBed -i - -c 7,8 -o mean,sum | awk '{OFS="\t"} {print $1, $2, $3,$4, $3-$2, $5/($3-$2)}' | sed 's/,/./g' > $OUTPUT/$SAMPLE".CG-intron"
intersectBed -b $C/$SAMPLE"_leaf.C-DMR.bed" -a $ANNO/"intron-wo-TE.merged.anno" -wao | awk '{if ($8>=99){print}}'|mergeBed -i - -c 7,8 -o mean,sum | awk '{OFS="\t"} {print $1, $2, $3,$4, $3-$2, $5/($3-$2)}' | sed 's/,/./g'  > $OUTPUT/$SAMPLE".C-intron"
intersectBed -b $CG/$SAMPLE"_leaf.CG-DMR.bed"  -a $ANNO/"intron-TE.merged.anno" -wao | awk '{if ($8>=99){print}}'| mergeBed -i - -c 7,8 -o mean,sum | awk '{OFS="\t"} {print $1, $2, $3,$4, $3-$2, $5/($3-$2)}' | sed 's/,/./g' > $OUTPUT/$SAMPLE".CG-intron-TE"
intersectBed -b $C/$SAMPLE"_leaf.C-DMR.bed" -a $ANNO/"intron-TE.merged.anno" -wao | awk '{if ($8>=99){print}}'|mergeBed -i - -c 7,8 -o mean,sum | awk '{OFS="\t"} {print $1, $2, $3,$4, $3-$2, $5/($3-$2)}' | sed 's/,/./g'  > $OUTPUT/$SAMPLE".C-intron-TE"

intersectBed -a $OUTPUT/$SAMPLE".CG-intron" -b $OUTPUT/$SAMPLE".C-intron" | awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}' > $OUTPUT/$SAMPLE"_PE-intron.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-intron" -b $OUTPUT/$SAMPLE"_PE-intron.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4="gbM"}'> $OUTPUT/$SAMPLE"_only-CG-intron.bed"
intersectBed -a $OUTPUT/$SAMPLE".C-intron" -b $OUTPUT/$SAMPLE"_PE-intron.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}'> $OUTPUT/$SAMPLE"_only-C-intron.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-intron-TE" -b $OUTPUT/$SAMPLE".C-intron-TE" | awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}' > $OUTPUT/$SAMPLE"_PE-intron-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-intron-TE" -b $OUTPUT/$SAMPLE"_PE-intron-TE.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4="gbM"}'> $OUTPUT/$SAMPLE"_only-CG-intron-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE".C-intron-TE" -b $OUTPUT/$SAMPLE"_PE-intron-TE.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4="teM"}'> $OUTPUT/$SAMPLE"_only-C-intron-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-intron" -b $OUTPUT/$SAMPLE".C-intron" | awk '{OFS="\t"}{print $1, $2, $3, $4}' > $OUTMETH/$SAMPLE"_PE-intron.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-intron" -b $OUTPUT/$SAMPLE"_PE-intron.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/$SAMPLE"_only-CG-intron.bed"
intersectBed -a $OUTPUT/$SAMPLE".C-intron" -b $OUTPUT/$SAMPLE"_PE-intron.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/$SAMPLE"_only-C-intron.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-intron-TE" -b $OUTPUT/$SAMPLE".C-intron-TE" | awk '{OFS="\t"}{print $1, $2, $3, $4}' > $OUTMETH/$SAMPLE"_PE-intron-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE".CG-intron-TE" -b $OUTPUT/$SAMPLE"_PE-intron-TE.bed" -v | awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/$SAMPLE"_only-CG-intron-TE.bed"
intersectBed -a $OUTPUT/$SAMPLE".C-intron-TE" -b $OUTPUT/$SAMPLE"_PE-intron-TE.bed" -v |  awk '{OFS="\t"}{print $1, $2, $3, $4}'> $OUTMETH/$SAMPLE"_only-C-intron-TE.bed"

cat $OUTPUT/$SAMPLE*"-gene.bed" | sort -k1,1n -k2,2n |awk '{OFS="\t"}{print}' > $WD/$SAMPLE"_epialleles-over-genes.bed"
cat $OUTPUT/$SAMPLE*"-gene-TE.bed" | sort -k1,1n -k2,2n | awk '{OFS="\t"}{print}' > $WD/$SAMPLE"_epialleles-over-genes-TE.bed"
cat $OUTPUT/$SAMPLE*"-intron.bed" | sort -k1,1n -k2,2n |awk '{OFS="\t"}{print}' > $WD/$SAMPLE"_epialleles-over-intron.bed"
cat $OUTPUT/$SAMPLE*"-exon.bed" | sort -k1,1n -k2,2n |awk '{OFS="\t"}{print}' > $WD/$SAMPLE"_epialleles-over-exon.bed"
cat $OUTPUT/$SAMPLE*"-intron-TE.bed" | sort -k1,1n -k2,2n |awk '{OFS="\t"}{print}' > $WD/$SAMPLE"_epialleles-over-intron-TE.bed"
cat $OUTPUT/$SAMPLE*"-exon-TE.bed" | sort -k1,1n -k2,2n |awk '{OFS="\t"}{print}' > $WD/$SAMPLE"_epialleles-over-exon-TE.bed"

cat $OUTMETH/$SAMPLE*"-gene.bed" | sort -k1,1n -k2,2n |awk '{OFS="\t"}{print}' > $OUTMETH/$SAMPLE"_epialleles-over-genes.bed"
cat $OUTMETH/$SAMPLE*"-gene-TE.bed" | sort -k1,1n -k2,2n | awk '{OFS="\t"}{print}' > $OUTMETH/$SAMPLE"_epialleles-over-genes-TE.bed"
cat $OUTMETH/$SAMPLE*"-intron.bed" | sort -k1,1n -k2,2n |awk '{OFS="\t"}{print}' > $OUTMETH/$SAMPLE"_epialleles-over-intron.bed"
cat $OUTMETH/$SAMPLE*"-exon.bed" | sort -k1,1n -k2,2n |awk '{OFS="\t"}{print}' > $OUTMETH/$SAMPLE"_epialleles-over-exon.bed"
cat $OUTMETH/$SAMPLE*"-intron-TE.bed" | sort -k1,1n -k2,2n |awk '{OFS="\t"}{print}' > $OUTMETH/$SAMPLE"_epialleles-over-intron-TE.bed"
cat $OUTMETH/$SAMPLE*"-exon-TE.bed" | sort -k1,1n -k2,2n |awk '{OFS="\t"}{print}' > $OUTMETH/$SAMPLE"_epialleles-over-exon-TE.bed"

