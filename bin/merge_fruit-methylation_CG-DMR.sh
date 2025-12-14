CHR="$( cat $1 )"
echo "------" $CHR
SAMPLE_DIR="/08_fruit-processing/08.1_DMR-classification/08.2_merge-DMRs/aa_get-ind-methylation-over-collapsed-loci"
CHR_SIZE="~/chr.size.bed"
OUT_DIR="08_fruit-processing/08.1_DMR-classification/08.2_merge-DMRs/ac_merge-methylation"
SAMPLE="$( ls ${SAMPLE_DIR}/*${CHR}_CG-DMR.methylation.bed | tr '\n' '\t' )"
LOCI_COLLAPSED="08_fruit-processing/08.1_DMR-classification/08.1_get-collapsed-loci"

#echo $SAMPLE
bedtools unionbedg -i $SAMPLE -empty -g ${CHR_SIZE} -filler NA |\
intersectBed -a - -b ${LOCI_COLLAPSED}/${CHR}_CG-DMR-loci_collapsed.bed > ${OUT_DIR}/${CHR}_CG-DMR.merged.methylation.bed

ls "${SAMPLE_DIR}"/*"${CHR}"*_CG_DMR.methylation.bed |\
    xargs -n 1 basename |\
    cut -d'_' -f1,2 > "${OUT_DIR}/00_${CHR}_CG-DMR_colNames.tsv"
