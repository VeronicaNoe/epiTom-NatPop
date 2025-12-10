#!/bin/bash
SAMPLE="$( cat $1 )"
CHR="$( cat $1 | cut -d'_' -f5 )"

SAMPLE_DIR="/mnt/disk2/vibanez/03_bismark_alignment/03.3_filtering/ac_filter/bb_output"
LOCI_DIR="/mnt/disk2/vibanez/08_fruit-processing/08.1_DMR-classification/08.1_get-collapsed-loci"
OUT_DIR="/mnt/disk2/vibanez/08_fruit-processing/08.1_DMR-classification/08.2_merge-DMRs/ab_weight-methylation/aa_data"

#extrar la metilacion en cada ventana definida con todas las muestras
for file in $(ls ${SAMPLE_DIR}/${SAMPLE}*$CHR*.filtered.bed); do
        base_name=$(basename "$file" .bed)
        cat ${file} | gawk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
        intersectBed -wb -a - -b ${LOCI_DIR}/$CHR"_C-DMR-loci_collapsed.bed" -nonamecheck |\
        gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum > ${OUT_DIR}/${base_name}.bed
done
