#!/bin/bash
READ_DIR="/mnt/disk2/vibanez/01_raw-data/01.0_biseq-reads"
SAMPLE="$( cat "$1" )"
OUT_DIR= $3

if [ "$2" == "ref" ]; then
	REF_DIR="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation"
else
	SAMPLE_PSEUDO_GENOME="$( cat "$1" | cut -d'_' -f1 )"
	REF_DIR="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/${SAMPLE_PSEUDO_GENOME}"
fi

bismark $REF_DIR \
        -1 ${READ_DIR}/${SAMPLE}"_1.f"*".gz" \
        -2 ${READ_DIR}/${SAMPLE}"_2.f"*".gz" \
        --bowtie2 \
        --phred33-quals \
        -p 8 \
        -N 1 \
        --un \
        --non_bs_mm \
        --ambiguous \
        -B ${SAMPLE} \
        -o ${OUT_DIR} \
        --temp_dir /mnt/disk1/vibanez/tmp
