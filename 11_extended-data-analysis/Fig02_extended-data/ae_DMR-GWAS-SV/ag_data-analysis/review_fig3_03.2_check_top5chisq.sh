#!/bin/bash
# Define the key file and output directory
KEY_FILE="/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results/03.1_SNPs_over_epiGenes_list.tsv"
OUTPUT_DIR="/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results/associated-SNP-over-epigenes"

# Preprocess the key file to create a regex pattern for zgrep
KEYS=$(paste -sd'|' "$KEY_FILE")

# Export the variables for parallel processing
export KEYS
export OUTPUT_DIR

# Use GNU parallel to process multiple files simultaneously
find /mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results/sig -name "*.ps.gz" | parallel -j80 '
  FILE={}
  OUTPUT_FILE="${OUTPUT_DIR}/$(basename "${FILE%.ps.gz}.snpsOverEpiGenes")"
  zgrep -E "^($KEYS)\b" "$FILE" > "$OUTPUT_FILE"
'

