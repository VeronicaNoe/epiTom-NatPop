#!/bin/bash
SAMPLE="$( cat $1 )"
IN_DIR="/mnt/disk2/vibanez/06_get-meth-vcf/aa_output"
OUT_DIR="/mnt/disk2/vibanez/06_get-meth-vcf/ab_epialleles"

# Find the input files
files=(${IN_DIR}/ch*_${SAMPLE}.vcf)
echo "Files to compress and index:"
printf "%s\n" "${files[@]}"

# Compress each file
parallel 'bgzip -c {} > {}.gz' ::: "${files[@]}"

# Index each compressed file
parallel 'tabix {}' ::: "${files[@]/%/.gz}"

# Concatenate with bcftools
bcftools concat "${files[@]/%/.gz}" -O z -o ${OUT_DIR}/${SAMPLE}.vcf.gz

# Clean up intermediate compressed files and indexes
#rm -f "${files[@]/%/.gz}" "${files[@]/%/.gz.tbi}"
