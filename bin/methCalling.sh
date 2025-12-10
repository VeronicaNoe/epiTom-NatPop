#!/bin/bash
SAMPLE="$( cat "$1" )"
echo $SAMPLE
OUT_DIR=$3
BASE_DIR=$(dirname "$(realpath "$0")")
IN_DIR=$(find "${BASE_DIR}/.." -type d \( -name "ab_deduplication" -o -name "03.2_deduplication" \) -print -quit)

if [ "$2" == "ref" ]; then
	REF_DIR="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation"
else
	SAMPLE_PSEUDO_GENOME="$( cat "$1" | cut -d'_' -f1 )"
	REF_DIR="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/${SAMPLE_PSEUDO_GENOME}"
fi

bismark_methylation_extractor \
	-p \
	--no_overlap \
	--comprehensive \
	--multicore 20 \
	--gzip \
	--gazillion \
	--CX \
	--buffer_size 60G \
	--cytosine_report \
	--CX \
      	--genome_folder ${REF_DIR} \
      	-o ${OUT_DIR} \
      	${IN_DIR}/$SAMPLE.sorted.deduplicated.bam

#coverage2cytosine \
#	--CX \
#	--genome_folder $REFDIR/$SAMPLE  \
#	-o $COVERAGE/$SAMPLE.with_cov \
#	$COVERAGE/$SAMPLE.*cov.gz
