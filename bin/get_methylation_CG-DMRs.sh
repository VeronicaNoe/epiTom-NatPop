SAMPLE="$( cat $1 | cut -d '_' -f1,2 )"
CHR="$( cat $1 | cut -d'_' -f3 )"

SAMPLE_DIR="/mnt/disk2/vibanez/03_bismark_alignment/03.3_filtering/ac_filter/bb_output"
LOCI_DIR="/mnt/disk2/vibanez/05_post-processing-DMRs/05.1_DMR-classification/05.1_get-collapsed-loci"
OUT_DIR="/mnt/disk2/vibanez/05_post-processing-DMRs/05.1_DMR-classification/05.2_merge-DMRs/aa_get-ind-methylation-over-collapsed-loci/ab_output"

#extrar la metilacion en cada ventana definida con todas las muestras
cat ${SAMPLE_DIR}/${SAMPLE}*_CG_${CHR}.filtered.bed | gawk '{OFS="\t"}{print $1,$2,$2,$4,$5}' |\
intersectBed -wb -a - -b $LOCI_DIR/$CHR"_CG-DMR-loci_collapsed.bed" -nonamecheck |\
gawk '{OFS="\t"}{print $1,$7,$8,$4,$5}' | mergeBed -c 4,5 -o sum,sum | gawk '{OFS="\t"}{print $1,$2,$3,$4/($4+$5)}' |\
sed 's/,/./g' > ${OUT_DIR}/${SAMPLE}_${CHR}_CG-DMR.methylation.bed

