#!/bin/bash
READ_DIR="01_raw-data/01.0_biseq-reads"
SAMPLE="$( cat "$1" )"
REF_DIR="03_biseq-processing/03.6_conversion-rate/aa_genome-preparation"
OUT_DIR="03_biseq-processing/03.6_conversion-rate/ab_mapping"

bismark $REF_DIR \
        -1 ${READ_DIR}/${SAMPLE}"_Biseq_1.f"*".gz" \
        -2 ${READ_DIR}/${SAMPLE}"_Biseq_2.f"*".gz" \
        --bowtie2 \
        --phred33-quals \
        -p 8 \
        -N 1 \
        --un \
        --non_bs_mm \
        --ambiguous \
        -B ${SAMPLE} \
        -o ${OUT_DIR} \
        --temp_dir tmp/
