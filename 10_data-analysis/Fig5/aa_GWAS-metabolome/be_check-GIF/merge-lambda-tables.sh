#!/bin/bash
GIF_DIR="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/be_check-GIF/cc_lambda-tables"
OUT_DIR="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/be_check-GIF"
# Set the working directory
cd ${GIF_DIR}
# Output file
lambda_file=${OUT_DIR}/"01.0_nonSig_by_lambda"

# Start with the header (copy from the first file)
head -n 1 $(ls *_lambda.tsv | head -n 1) > "$output"

# Then append all content (skip header lines)
for file in *_lambda.tsv; do
    tail -n +2 "$file" >> ${lambda_file}
echo "✅ All lambda tables merged into: $lambda_file"
done
