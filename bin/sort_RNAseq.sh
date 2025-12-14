#!/bin/bash
SAMPLE="$( cat "$1" )"
IN_DIR="07_rnaseq-processing/07.1_mapping"
OUT_DIR="07_rnaseq-processing/07.2_sorting"

samtools sort -@20 -T /tmp $IN_DIR/$SAMPLE.bam -o $OUT_DIR/$SAMPLE.sorted.bam -O BAM
samtools index -@20 $OUT_DIR/$SAMPLE.sorted.bam
