#!/bin/bash
# conda activate tools
SAMPLE="$( echo "$1" | sed 's/.sorted.bam//g' )"
KODIR="/mnt/disk2/vibanez/02_methylkit/am_KO-RNAseq/counts"
#OUTDIR="/mnt/disk2/vibanez/tomato_RNAseq/03_counts"
ANNO="/mnt/disk2/vibanez/tomato_RNAseq/gtf_files/sol2.4_nonRef.gff3"

featureCounts -a "$ANNO" \
	-o $KODIR/$SAMPLE".counts.tsv" \
	--largestOverlap \
	-t exon \
	-g Parent \
	-F "GFF3" \
	$SAMPLE".sorted.bam"

