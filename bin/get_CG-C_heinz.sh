#!/bin/bash
CHR="$( cat $1 )"
echo $CHR
SAMPLE_PATH="/mnt/disk2/vibanez/03_biseq-processing/03.3_filtering/ac_filter/bc_heinz"
TMP="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/tmp-files"
LOCI_COLLAPSED="/mnt/disk2/vibanez/05_post-processing-DMRs/05.1_DMR-classification/05.1_get-collapsed-loci"
OUT_DIR="/mnt/disk2/vibanez/05_post-processing-DMRs/05.1_DMR-classification/05.2_merge-DMRs/aa_get-ind-methylation-over-collapsed-loci/ab_output"

# intersect windows and samples
echo "------ get tmp C-DMR"
cat $SAMPLE_PATH/"TS-253_leaf_Biseq_C"*$CHR".filtered.bed" | awk '{OFS="\t"}{print $1, $2, $2+99, $4,$5}' |\
intersectBed -a - -b $LOCI_COLLAPSED/$CHR"_C-DMR-loci_collapsed.bed" -wb |\
awk '{OFS="\t"}{print $6, $7, $8, $4,$5}'| sortBed -i - | mergeBed -i - -c 4,5 -o sum,sum |\
awk '{OFS="\t"}{print $1, $2, $3, $4/($4+$5)}'  > $TMP/$CHR"_TS-253_leaf_C-DMR.tmp"

echo "------ get CG-DMR"
cat $SAMPLE_PATH/"TS-253_leaf_Biseq_CG_"$CHR".filtered.bed" | awk '{OFS="\t"}{print $1, $2, $2+99, $4,$5}' |\
intersectBed -a - -b $LOCI_COLLAPSED/$CHR"_CG-DMR-loci_collapsed.bed" -wb |\
awk '{OFS="\t"}{print $6, $7, $8, $4,$5}' | sortBed -i - | mergeBed -i - -c 4,5 -o sum,sum |\
bedtools subtract -a - -b  $TMP/$CHR"_TS-253_leaf_C-DMR.tmp" |\
awk '{OFS="\t"}{print $1, $2, $3, $4/($4+$5)}' > $OUT_DIR/$CHR"_TS-253_leaf.CG-DMR.bed"

echo "------ get C-DMR"
intersectBed -a $TMP/$CHR"_TS-253_leaf_C-DMR.tmp" -b $OUT_DIR/$CHR"_TS-253_leaf.CG-DMR.bed" -v > $OUT_DIR/$CHR"_TS-253_leaf.C-DMR.bed"
