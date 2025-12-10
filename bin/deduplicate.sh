#!/bin/bash
SAMPLE="$( cat "$1" )"
echo $SAMPLE
MAP_DIR="/mnt/disk2/vibanez/03_biseq-processing/03.1_mapping"
OUT_DIR=$2
#"/mnt/disk2/vibanez/03_biseq-processing/03.2_deduplication"
# sort
samtools sort -n -@30 -T /mnt/disk1/vibanez/tmp/ ${MAP_DIR}/${SAMPLE}*".bam" -o ${OUT_DIR}/$SAMPLE.sorted.bam -O BAM
# deduplicate
deduplicate_bismark -p --bam ${OUT_DIR}/$SAMPLE.sorted.bam -o ${OUT_DIR}/$SAMPLE.sorted
