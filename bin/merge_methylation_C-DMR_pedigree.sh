#!/bin/bash
GEN="$( cat $1 )"
SAMPLE_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ab_pedigree/aa_get-ind-methylation-over-collapsed-loci"
LOCI="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci/ab_pedigree"
CHR_SIZE="/home/vibanez/bin/chr.size.bed"
OUT_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ab_pedigree/ab_merge-methylation"
SAMPLE="$( ls ${SAMPLE_DIR}/*"_C-DMR_"${GEN}".methylation.bed" | tr '\n' '\t' )"

echo $SAMPLE
bedtools unionbedg -i ${SAMPLE} -empty -g ${CHR_SIZE} -filler NA |\
intersectBed -a - -b ${LOCI}/"C-DMR-loci_collapsed.bed" > ${OUT_DIR}/"C-DMR_"${GEN}".merged.methylation.bed"

ls "${SAMPLE_DIR}"/*"_C-DMR_all.methylation.bed" |\
    xargs -n 1 basename |\
    cut -d'_' -f1 > "${OUT_DIR}/C-DMR_"${GEN}".colNames.tsv"
