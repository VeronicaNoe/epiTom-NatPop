#!/bin/bash
SAMPLE="$( cat $1 )"
echo $SAMPLE
CHR="$( cat $1 | cut -d'_' -f5 )"
echo $CHR
SAMPLEPATH="/mnt/disk2/vibanez/01_bismark_analysis/04_processingFiles/ac_filter/bb_output"
C_LOCI_DIR="/mnt/disk2/vibanez/02_methylkit/aj_common/loci"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/aj_common/ac_merge_C-DMR/aa_data"

cat $SAMPLEPATH/$SAMPLE.filtered.bed | awk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
sortBed -i - | intersectBed -wb -a - -b $C_LOCI_DIR/$CHR"_C-DMR-loci_collapsed.bed" -nonamecheck | awk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | \
mergeBed -c 4,5 -o sum,sum | awk '{OFS="\t"}{print $1,$2,$3,$4,$5}'  > $OUTPATH/$SAMPLE.methylation.bed
