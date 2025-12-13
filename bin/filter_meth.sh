#!/bin/bash
SAMPLE="$( cat "$1" )"
SAMPLE_DIR="03_biseq-processing/03.4_filtering/aa_ctxt-split"
OUT_DIR="03_biseq-processing/03.4_filtering/ab_chr-split/bb_data"
#<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
echo $SAMPLE
cat ${SAMPLE_DIR}/${SAMPLE}.bed | gawk '{if ($4 + $5 >= 5) {print}}' > $OUT_DIR/$SAMPLE.prefiltered.bed
