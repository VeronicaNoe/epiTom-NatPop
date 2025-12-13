#!/bin/bash
CHR="$( cat $1 )"
SAMPLE_DIR="05_DMR-processing/05.1_DMR-classification/05.0_get-DMR/ab_output"
OUT_DIR="05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci"
TMP="tmp"
# creating a small data set
for file in $(ls ${SAMPLE_DIR}/$CHR*.C-DMR.bed); do
        base_name=$(basename "$file" .bed)
        cat $file | gawk '{OFS="\t"}{print $1,$2,$3,$6/$5}' | sed 's/,/./g' > $TMP/${base_name}.tmp
done
#count the overlaps
cat $TMP/$CHR*.C-DMR.tmp | sort -k 1.4n -k 2n | gawk '{OFS="\t"}{print $1,$2,$3,$4}'|\
mergeBed -c 4,4 -o count,distinct | gawk '{OFS="\t"}{print $1,$2,$3,$4,$5}'| uniq > ${OUT_DIR}/$CHR"_C-DMR-loci_collapsed.bed"

#loop over all the files to have exactly the same coordinates for all, needed to run multiinter
for FILE in ${SAMPLE_DIR}/$CHR*.C-DMR.bed; do
    base_name=$(basename "$FILE" .bed)
    awk '{OFS="\t"}{print $1,$2,$3}' "$FILE" |\
    intersectBed -wa -a ${OUT_DIR}/$CHR"_C-DMR-loci_collapsed.bed" -b - -nonamecheck | uniq > ${OUT_DIR}/ab_output/${base_name}.loci.bed
done
#rm $TMP/$CHR*.C-DMR.tmp
