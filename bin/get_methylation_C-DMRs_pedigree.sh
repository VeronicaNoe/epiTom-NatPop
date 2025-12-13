#!/bin/bash
SAMPLE="$( cat $1 | cut -d '_' -f1 )"
SAMPLE_DIR="03_biseq-processing/03.4_filtering/ac_filter"
LOCI_DIR="05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci/ab_pedigree"
OUT_DIR="05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ab_pedigree/aa_get-ind-methylation-over-collapsed-loci"

#extrar la metilacion en cada ventana definida con todas las muestras
cat ${SAMPLE_DIR}/${SAMPLE}_CH*.filtered.bed | gawk '{OFS="\t"}{print $1,$2-1, $2,$4, $5}' |\
intersectBed -wb -a - -b $LOCI_DIR/"C-DMR-loci_collapsed.bed" -nonamecheck |\
sortBed -i - |\
gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
sed 's/,/./g' > ${OUT_DIR}/${SAMPLE}_C-DMR_all.methylation.bed


#cat ${SAMPLE_DIR}/${SAMPLE}*_*.filtered.bed | gawk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
#intersectBed -wb -a - -b $LOCI_DIR/"C-DMR_interGen.positions" -nonamecheck |\
#gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | sortBed -i - | mergeBed -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
#sed 's/,/./g' > ${OUT_DIR}/${SAMPLE}_C-DMR_interGen.methylation.bed

#cat ${SAMPLE_DIR}/${SAMPLE}*_*.filtered.bed | gawk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
#intersectBed -wb -a - -b $LOCI_DIR/"C-DMR_intraG0.positions" -nonamecheck |\
#gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | sortBed -i - | mergeBed -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
#sed 's/,/./g' > ${OUT_DIR}/${SAMPLE}_C-DMR_intraG0.methylation.bed

#cat ${SAMPLE_DIR}/${SAMPLE}*_*.filtered.bed | gawk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
#intersectBed -wb -a - -b $LOCI_DIR/"C-DMR_intraG5.positions" -nonamecheck |\
#gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | sortBed -i - | mergeBed -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
#sed 's/,/./g' > ${OUT_DIR}/${SAMPLE}_C-DMR_intraG5.methylation.bed

