#!/bin/bash
SAMPLE="$( cat "$1" )"
IN_DIR="07_rnaseq-processing/07.2_sorting"
OUT_DIR="07_rnaseq-processing/07.3_counts"
ANNO_DIR="07_rnaseq-processing/07.0_ref/gtf_files/sol2.4_nonRef.gff3"

featureCounts -p -a $ANNO_DIR \
	-o $IN_DIR/$SAMPLE".counts.tsv" \
	--largestOverlap \
	-t exon \
	-g Parent \
	-F "GFF3" \
	$SAMPLE".sorted.bam"

cat ${IN_DIR}/${SAMPLE}".counts.tsv" | grep -v "#" |awk '{OFS="\t"}{print $1,$7}' | sed 's/.sorted.bam//g' |\
  sed 's/mRNA://g' > ${IN_DIR}/${SAMPLE}"_edited.tsv"
