#!/bin/bash
OUT_DIR="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/be_check-GIF"
sig_dir="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results/sig"
nonsig_dir="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results/nonSig-GIF"
lambda_file=${OUT_DIR}/"01.0_nonSig_by_lambda"

# Create nonSig-GIF folder if it doesn't exist
mkdir -p "$nonsig_dir"

# Read the file line-by-line, skipping header
tail -n +2 "$lambda_file" | while read -r sample metabolite molmarker kinship index number_QLTs lambda
do
    # Check if lambda is outside [0.975, 1.025]
    result=$(awk -v l="$lambda" 'BEGIN {if (l <= 0.95 || l >= 1.05) print "move"; else print "ok"}')

    if [ "$result" == "move" ]; then
        echo "Moving: $sample.ps.gz (lambda = $lambda)"
        mv $sig_dir/${sample}.* $nonsig_dir/
    else
        echo "OK: $sample.ps.gz (lambda = $lambda)"
    fi
done
