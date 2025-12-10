#!/bin/bash
READ_DIR="/mnt/disk2/vibanez/01_raw_data/01.0_biseq-reads"
REF_DIR="/mnt/disk2/vibanez/02_pseudo-pan-genomes/02.1_make-pseudogenomes"
OUT_DIR="/mnt/disk2/vibanez/03_bismark_alignment/03.1_mapping"

SAMPLE="$( cat "$1" )"
echo $SAMPLE

bismark $REF_DIR/$SAMPLE \
        ${READ_DIR}/$SAMPLE".f"*".gz" \
        --bowtie2 \
        --phred33-quals \
        -p 8 \
        -N 1 \
        --un \
        --non_bs_mm \
        --ambiguous \
        -B $SAMPLE \
        -o $OUT_DIR \
        --temp_dir /mnt/disk1/vibanez/tmp
