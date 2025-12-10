#!/bin/bash
READ_DIR="/mnt/disk2/vibanez/reads"
SAMPLE="$( cat "$1" )"
REF_DIR="/mnt/disk6/vibanez/conversion-rate/aa_genome-preparation"
OUT_DIR="/mnt/disk6/vibanez/conversion-rate/ab_mapping"

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
        --temp_dir /mnt/disk6/vibanez/tmp/
