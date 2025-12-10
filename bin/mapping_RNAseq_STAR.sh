#!/bin/bash
SAMPLE="$(cat "$1")"
OUT_DIR="$2"
IN_DIR="/mnt/disk2/vibanez/01_raw-data/01.1_rnaseq-reads"
REF_DIR="/mnt/disk2/vibanez/07_rnaseq-processing/07.0_ref"

# Detect single-end or paired-end based on available files
if ls $IN_DIR/${SAMPLE}*"_se."*".gz" 1> /dev/null 2>&1; then
    STAR --runThreadN 60 \
        --runMode alignReads \
        --genomeDir $REF_DIR/TS-253 \
        --readFilesIn $IN_DIR/${SAMPLE}*"_se"*".gz" \
        --readFilesCommand zcat \
        --quantMode GeneCounts \
        --outFilterMultimapNmax 100 \
        --winAnchorMultimapNmax 100 \
        --outFilterMismatchNmax 2 \
        --outMultimapperOrder Random \
        --outFileNamePrefix $SAMPLE \
        --outSAMmultNmax 1 \
        --outSAMattributes All \
        --twopass1readsN -1 \
        --twopassMode Basic \
        --outSAMtype BAM Unsorted \
        --outStd BAM_Unsorted > $OUT_DIR/$SAMPLE.bam
else
    STAR --runThreadN 60 \
        --runMode alignReads \
        --genomeDir $REF_DIR/TS-253 \
        --readFilesIn $IN_DIR/${SAMPLE}*"_1."*".gz" $IN_DIR/$SAMPLE*"_2."*".gz" \
        --readFilesCommand zcat \
        --quantMode GeneCounts \
        --outFilterMultimapNmax 100 \
        --winAnchorMultimapNmax 100 \
        --outFilterMismatchNmax 2 \
        --outMultimapperOrder Random \
        --outFileNamePrefix $SAMPLE \
        --outSAMmultNmax 1 \
        --outSAMattributes All \
        --twopass1readsN -1 \
        --twopassMode Basic \
        --outSAMtype BAM Unsorted \
        --outStd BAM_Unsorted > $OUT_DIR/$SAMPLE.bam
fi

# Move the report
mv $OUT_DIR/${SAMPLE}_report.txt $OUT_DIR/*_reports/
