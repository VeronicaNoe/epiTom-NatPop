#!/bin/bash
DIR="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.1_get-collapsed-loci/ab_pedigree"

cat ${DIR}/*_CG-DMR-loci_collapsed.bed | sortBed -i - | mergeBed -i - > ${DIR}/"CG-DMR_all.positions"
cat ${DIR}/*_C-DMR-loci_collapsed.bed | sortBed -i - | mergeBed -i - > ${DIR}/"C-DMR_all.positions"
