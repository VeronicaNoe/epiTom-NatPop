#!/bin/bash
IN_DIR="/mnt/disk2/vibanez/03_biseq-processing/03.4_filtering/ac_filter/5x"
OUT_DIR="/mnt/disk2/vibanez/03_biseq-processing/03.5_igv-files"
SAMPLE="$( cat "$1" )"
echo $SAMPLE
# --- Commands ---
#echo ls $SAMPLE*bed
echo '---------- Merging chr sample: ' ${SAMPLE}
cat $IN_DIR/${SAMPLE}*.filtered.bed | awk '{OFS="\t"}{print $1,$2-1,$2,$4,$5}' |\
sortBed -i - | mergeBed -i - -o sum,sum -c 4,5 |\
awk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' | grep 'SL2.50' > $OUT_DIR/${SAMPLE}.bed
echo '---------- Formating chr sample: ' ${SAMPLE}
~/bin/bedGraphToBigWig ${OUT_DIR}/${SAMPLE}.bed ~/bin/chr.size.bed $OUT_DIR/$SAMPLE.bw
#rm ${OUT_DIR}/${SAMPLE}.bed
