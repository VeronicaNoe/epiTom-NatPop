#!/bin/bash
SAMPLE="$( cat "$1" )"
echo $SAMPLE
REF_DIR="/mnt/disk6/vibanez/conversion-rate/aa_genome-preparation"
IN_DIR="/mnt/disk6/vibanez/conversion-rate/ab_mapping"
OUT_DIR="/mnt/disk6/vibanez/conversion-rate/ac_methylation-calling"
bismark_methylation_extractor \
	-p \
	--no_overlap \
	--comprehensive \
	--multicore 20 \
	--gzip \
	--gazillion \
	--CX \
	--buffer_size 100G \
	--CX \
	--cytosine_report \
      	--genome_folder ${REF_DIR} \
      	-o ${OUT_DIR} \
      	${IN_DIR}/$SAMPLE.bam

zcat ${OUT_DIR}/${SAMPLE}.CX_report.txt.gz | grep 'CHH' | \
awk '{methylated += $4; unmethylated += $5} END {print "'$SAMPLE'," ,1 - (methylated/(methylated+unmethylated))}' \
> ${OUT_DIR}/${SAMPLE}_conversion_rate.txt
