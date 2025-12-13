#!/bin/bash
CHR="$( cat $1 )"
SAMPLE_DIR="05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/aa_get-ind-methylation-over-collapsed-loci"
CHR_SIZE="~/bin/chr.size.bed"
OUT_DIR="05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation"
SAMPLE="$( ls ${SAMPLE_DIR}/$CHR*.weighted.methylation.bed | tr '\n' '\t' )"
LOCI_COLLAPSED="05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci/aa_natural-accessions"
#echo $SAMPLE

bedtools unionbedg -i $SAMPLE -empty -g $CHR_SIZE -filler NA |\
intersectBed -a - -b ${LOCI_COLLAPSED}/$CHR"_C-DMR-loci_collapsed.bed" > ${OUT_DIR}/$CHR"_C-DMR.merged.methylation.bed"

ls "${SAMPLE_DIR}"/*"${CHR}"*_CG_DMR.methylation.bed |\
    xargs -n 1 basename |\
    cut -d'_' -f1,2 > "${OUT_DIR}/00_${CHR}_CG-DMR_colNames.tsv"

