#!/bin/bash
OUT_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci/ab_pedigree"
TMP="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.0_get-DMR/ab_predigree/ctxt-dmr"

cat $TMP/*CH*DMR_merged.bed | sortBed -i - | gawk '{OFS="\t"}{print $1,$2,$3}'|\
mergeBed -i - > ${OUT_DIR}/"C-DMR-loci_collapsed.bed"

