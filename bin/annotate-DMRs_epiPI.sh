#!/bin/bash
GROUP="$( cat $1 | cut -d'_' -f1 )"
ANNO="$( cat $1 | cut -d'_' -f2 )"
echo $ANNO
INDIR="/mnt/data6/vibanez/mVCF/ac_calculate-ePI/ba_data"
WD_DIR="/mnt/data6/vibanez/mVCF/ac_calculate-ePI/data"
OUTDIR="/mnt/data6/vibanez/mVCF/ac_calculate-ePI/bb_output"
ANNODIR="/mnt/data6/vibanez/SNPs/vcfiles/ad_annotation/data"
DMR_BED="/mnt/data6/vibanez/mVCF/ac_calculate-ePI/data"
echo "CG-DMR"
## this does not have sense
## use the annotate-DMRs
#intersectBed -a $DMR_BED/CG-DMR.bed -b $ANNODIR/$ANNO.anno | sortBed -i -  | mergeBed -i - > $WD_DIR/"CG-DMR_"$ANNO.bed
## new version
awk 'NR>1' $INDIR/$GROUP"_CG-DMR_general_allAcc.recode.windowed.pi" | sortBed -i - |intersectBed -a - -b $WD_DIR/"CG-DMR_"$ANNO.bed -wa > $OUTDIR/"CG-DMR_"$GROUP"_"$ANNO"_allAcc.windowed.pi"

echo "C-DMR"
#intersectBed -a $DMR_BED/C-DMR.bed -b $ANNODIR/$ANNO.anno  | sortBed -i - | mergeBed -i - > $WD_DIR/"C-DMR_"$ANNO.bed
awk 'NR>1' $INDIR/$GROUP"_C-DMR_general_allAcc.recode.windowed.pi" | intersectBed -a - -b $WD_DIR/"C-DMR_"$ANNO.bed -wa > $OUTDIR/"C-DMR_"$GROUP"_"$ANNO"_allAcc.windowed.pi"

