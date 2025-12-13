#!/bin/bash
SAMPLE="$( cat $1 | cut -d '_' -f1,2 )"
CHR="$( cat $1 | cut -d'_' -f3 )"

SAMPLE_DIR="03_biseq-processing/03.3_filtering/ac_filter"
LOCI_DIR="05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci"
OUT_DIR="05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/aa_get-ind-methylation-over-collapsed-loci"

#extrar la metilacion en cada ventana definida con todas las muestras
for file in $(ls ${SAMPLE_DIR}/${SAMPLE}*$CHR*.filtered.bed); do
        base_name=$(basename "$file" .bed)
	cat ${file} | gawk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
	intersectBed -wb -a - -b $LOCI_DIR/$CHR"_C-DMR-loci_collapsed.bed" -nonamecheck |\
	gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum > ${OUT_DIR}/${base_name}.bed
done
