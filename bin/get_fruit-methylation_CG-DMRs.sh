SAMPLE="$( cat $1 )"
CHR="$( cat $1 | cut -d'_' -f5 )"

SAMPLE_DIR="03_biseq-processing/03.4_filtering/ac_filter"
LOCI_DIR="08_fruit-processing/08.1_DMR-classification/08.1_get-collapsed-loci"
OUT_DIR="08_fruit-processing/08.1_DMR-classification/08.2_merge-DMRs/aa_get-ind-methylation-over-collapsed-loci"

cat ${SAMPLE_DIR}/${SAMPLE}*_CG_${CHR}.filtered.bed | gawk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
intersectBed -wb -a - -b ${LOCI_DIR}/${CHR}"_CG-DMR-loci_collapsed.bed" -nonamecheck |\
awk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum | awk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
sed 's/,/./g' > ${OUT_DIR}/${SAMPLE}_${CHR}"_CG-DMR.methylation.bed"

