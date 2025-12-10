#!/bin/bash
LEAF_CG="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/ae_loci_CG-DMRs/tmp"
LEAF_C="/mnt/disk2/vibanez/02_methylkit/ab_preprocessing/af_loci_C-DMRs/tmp"
FRUIT_CG="/mnt/disk2/vibanez/02_methylkit/FRUITs/ad_loci_CG-DMRs/tmp"
FRUIT_C="/mnt/disk2/vibanez/02_methylkit/FRUITs/ae_loci_C-DMRs/tmp"

OUTPATH="/mnt/disk2/vibanez/02_methylkit/FRUITs/ai_tissue-merge/loci"

echo "=========  CG  ==========="
echo "--- leaf"
cat $LEAF_CG/CG-DMR-loci_collapsed.bed >> $OUTPATH/CG-DMR-loci_collapsed.tmp
echo "--- fruit"
cat $FRUIT_CG/CG-DMR-loci_collapsed.bed >> $OUTPATH/CG-DMR-loci_collapsed.tmp
echo "--- sorting CG"
sortBed -i $OUTPATH/CG-DMR-loci_collapsed.tmp | mergeBed -i - > $OUTPATH/CG-DMR-loci_collapsed.bed

echo "=========  C  ============"
echo "--- leaf"
cat $LEAF_C/C-DMR-loci_collapsed.bed >> $OUTPATH/C-DMR-loci_collapsed.tmp
echo "--- fruit"
cat $FRUIT_C/C-DMR-loci_collapsed.bed >> $OUTPATH/C-DMR-loci_collapsed.tmp
echo "--- sorting"
sortBed -i $OUTPATH/C-DMR-loci_collapsed.tmp | mergeBed -i - > $OUTPATH/C-DMR-loci_collapsed.bed

rm $OUTPATH/*.tmp
