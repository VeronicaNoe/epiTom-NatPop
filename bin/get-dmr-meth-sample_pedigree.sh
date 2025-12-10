#!/bin/bash
SAMPLE="$( cat $1 | cut -d '_' -f1 )"

SAMPLE_DIR="/mnt/disk2/vibanez/03_biseq-processing/03.4_filtering/ac_filter"
DMR_DIR="/mnt/disk2/vibanez/04_methylome-comparison/ab_pedigree"
TMP_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.0_get-DMR/ab_predigree/ctxt-dmr"

cat ${DMR_DIR}/*${SAMPLE}*CG_DMR.methylation.bed | gawk '{OFS="\t"}{print $1,$2,$3}' |\
sortBed -i - > ${TMP_DIR}/${SAMPLE}_CG-DMR_merged.bed
#extrar la metilacion en cada ventana definida con todas las muestras
cat ${SAMPLE_DIR}/${SAMPLE}_CG_*.filtered.bed | gawk '{OFS="\t"}{print $1,$2-1,$2,$4,$5}' |\
intersectBed -wb -a - -b ${TMP_DIR}/${SAMPLE}_"CG-DMR_merged.bed" -nonamecheck |\
sortBed -i - |\
gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1,$2,$3,$4,$5}' |\
sed 's/,/./g' > ${TMP_DIR}/${SAMPLE}_CG-DMR.methylation.bed

cat ${DMR_DIR}/*${SAMPLE}*CHG_DMR.methylation.bed | gawk '{OFS="\t"}{print $1,$2,$3}' |\
sortBed -i - > ${TMP_DIR}/${SAMPLE}_CHG-DMR_merged.bed
#extrar la metilacion en cada ventana definida con todas las muestras
cat ${SAMPLE_DIR}/${SAMPLE}_CHG_*.filtered.bed | gawk '{OFS="\t"}{print $1,$2-1,$2,$4,$5}' |\
intersectBed -wb -a - -b ${TMP_DIR}/${SAMPLE}_"CHG-DMR_merged.bed" -nonamecheck |\
sortBed -i - |\
gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1,$2,$3,$4,$5}' |\
sed 's/,/./g' > ${TMP_DIR}/${SAMPLE}_CHG-DMR.methylation.bed

cat ${DMR_DIR}/*${SAMPLE}*CHH_DMR.methylation.bed | gawk '{OFS="\t"}{print $1,$2,$3}' |\
sortBed -i - > ${TMP_DIR}/${SAMPLE}_CHH-DMR_merged.bed
#extrar la metilacion en cada ventana definida con todas las muestras
cat ${SAMPLE_DIR}/${SAMPLE}_CHH_*.filtered.bed | gawk '{OFS="\t"}{print $1,$2-1,$2,$4,$5}' |\
intersectBed -wb -a - -b ${TMP_DIR}/${SAMPLE}_"CHH-DMR_merged.bed" -nonamecheck |\
sortBed -i - |\
gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1,$2,$3,$4,$5}' |\
sed 's/,/./g' > ${TMP_DIR}/${SAMPLE}_CHH-DMR.methylation.bed
