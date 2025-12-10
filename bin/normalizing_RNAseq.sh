#!/bin/bash
SAMPLE="$( cat $1 )"
IN_DIR="/mnt/disk2/vibanez/07_rnaseq-processing/07.2_sorting"
OUT_DIR="/mnt/disk2/vibanez/07_rnaseq-processing/07.3_counts"

bamCoverage --numberOfProcessors 30 \
	--binSize 50 \
	--effectiveGenomeSize 1070311436 \
	--normalizeUsing RPGC \
	--exactScaling \
	--outFileFormat bedgraph \
	--bam $IN_DIR/$SAMPLE".sorted.bam" \
	-o $OUT_DIR/$SAMPLE"_normalized"
