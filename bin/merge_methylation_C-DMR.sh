#!/bin/bash
CHR="$( cat $1 )"
SAMPLE_DIR="/mnt/disk2/vibanez/05_post-processing-DMRs/05.1_DMR-classification/05.2_merge-DMRs/ab_weight-methylation/bc_output"
CHR_SIZE="~/bin/chr.size.bed"
OUT_DIR="/mnt/disk2/vibanez/05_post-processing-DMRs/05.1_DMR-classification/05.2_merge-DMRs/ac_merge-methylation/bb_output"
SAMPLE="$( ls ${SAMPLE_DIR}/$CHR*.weighted.methylation.bed | tr '\n' '\t' )"
LOCI_COLLAPSED="/mnt/disk2/vibanez/05_post-processing-DMRs/05.1_DMR-classification/05.1_get-collapsed-loci"
#echo $SAMPLE

bedtools unionbedg -i $SAMPLE -empty -g $CHR_SIZE -filler NA |\
intersectBed -a - -b ${LOCI_COLLAPSED}/$CHR"_C-DMR-loci_collapsed.bed" > ${OUT_DIR}/$CHR"_C-DMR.merged.methylation.bed"

ls "${SAMPLE_DIR}"/*"${CHR}"*_CG_DMR.methylation.bed |\
    xargs -n 1 basename |\
    cut -d'_' -f1,2 > "${OUT_DIR}/00_${CHR}_CG-DMR_colNames.tsv"

