#!/bin/bash
CHR="$( cat $1 )"
SAMPLE_DIR="/mnt/disk2/vibanez/08_fruit-processing/08.1_DMR-classification/08.0_get-DMR/bb_output"
OUT_DIR="/mnt/disk2/vibanez/08_fruit-processing/08.1_DMR-classification/08.1_get-collapsed-loci"
TMP="/mnt/disk2/vibanez/08_fruit-processing/08.1_DMR-classification/tmp-files"
# creating a small data set with all the windows
for file in $(ls ${SAMPLE_DIR}/$CHR*.bed); do
	base_name=(basename "$file" .bed)
	cat $file | gawk '{OFS="\t"}{print $1,$2,$3,$9/$8}' | sed 's/,/./g' > $TMP/${base_name}.tmp
done
#count the overlaps by each windows
cat $TMP/$CHR*.CG-DMR.tmp| sort -k 1.4n -k 2n | gawk '{OFS="\t"}{print $1,$2,$3,$4}' |\
mergeBed -c 4,4 -o count,distinct | gawk '{OFS="\t"}{print $1,$2,$3,$4,$5}'| uniq > ${OUT_DIR}/$CHR"_CG-DMR-loci_collapsed.bed"

#loop over all the files to have exactly the same coordinates for all, needed to run multiinter
for FILE in $(ls ${SAMPLE_DIR}/$CHR*.CG-DMR.bed ); do
	base_name=$(basename "$FILE" .bed)
	awk '{OFS="\t"}{print $1,$2,$3}' "$FILE" |\
	intersectBed -wa -a ${OUT_DIR}/$CHR"_CG-DMR-loci_collapsed.bed" -b - -nonamecheck | uniq > ${OUT_DIR}/ab_output/${base_name}.loci.bed
done
rm $TMP/$CHR*.CG-DMR.tmp
