#!/bin/bash
SAMPLE="$( cat "$1" )"
ROOT="/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/ac_get-DMR-ONT"
#SAMPLE_DIR="${ROOT}/02_ctxt-split"
SAMPLE_DIR="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/ae_ctxt-split"
#OUT_DIR="${ROOT}/03_chr-split/ab_data"
OUT_DIR="/mnt/disk2/vibanez/otherAnalysis/19_pedigree-denovo-assembly/af_chr-split"
mkdir -p ${OUT_DIR}
#<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
# we will filter in mkit
echo $SAMPLE
cat ${SAMPLE_DIR}/${SAMPLE}.bed | gawk '{if ($4 + $5 >= 1) {print}}' > $OUT_DIR/$SAMPLE.prefiltered.bed
