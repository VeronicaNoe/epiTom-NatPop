#!/bin/bash
GROUP="$( cat $1 | cut -d'_' -f1 )"
ANNO="$( cat $1 | cut -d'_' -f2 )"
echo $ANNO
#INDIR="/mnt/new_disk2/data6/vibanez/SNPs/vcfiles/aa_split-vcf-group/bb_output"
INDIR="/mnt/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bb_general"
WD_DIR="/mnt/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bc_annotated/data"
OUTDIR="/mnt/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bc_annotated/bb_output"
ANNODIR="/mnt/data6/vibanez/SNPs/vcfiles/ad_annotation/data"
SNP_BED="/mnt/data6/vibanez/SNPs/vcfiles/ab_calculate-PI/bc_annotated/data/SNPs.bed"

#intersectBed -a $SNP_BED -b $ANNODIR/$ANNO.anno  | mergeBed -i - > $WD_DIR/"SNPs_"$ANNO.bed
awk 'NR>1' $INDIR/$GROUP"_allAcc_SNPs.windowed.pi" | intersectBed -a - -b $WD_DIR/"SNPs_"$ANNO.bed -wa > $OUTDIR/"SNPs_"$GROUP"_"$ANNO"_allAcc.windowed.pi"
awk 'NR>1' $INDIR/$GROUP"_20acc_SNPs.windowed.pi" | intersectBed -a - -b $WD_DIR/"SNPs_"$ANNO.bed -wa > $OUTDIR/"SNPs_"$GROUP"_"$ANNO"_20acc.windowed.pi"

