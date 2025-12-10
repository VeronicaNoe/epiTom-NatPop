import os
import pandas as pd
import numpy as np
from scipy.stats import chi2

# Define directories
data_dir = "/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results/sig"
out_dir = "/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/be_check-GIF"

# List to store samples failing the lambda check
remove_by_lambda = []

# Ensure input directory exists
if not os.path.exists(data_dir):
    raise FileNotFoundError(f"Input directory does not exist: {data_dir}")

# Get list of input files
input_files = [f for f in os.listdir(data_dir) if f.endswith(".ps.gz")]

if not input_files:
    print(f"No .ps.gz files found in directory: {data_dir}")
    exit()

for file_name in input_files:
    file_path = os.path.join(data_dir, file_name)
    
    # Ensure file exists (redundant with os.listdir but safe for manual edits)
    if not os.path.isfile(file_path):
        print(f"File not found: {file_path}")
        continue

    try:
        # Read the data
        df = pd.read_csv(file_path, sep='\t', na_values=["NA"], compression='gzip', header=None)

        # Check if the file has the expected number of columns
        if len(df.columns) < 4:
            print(f"File {file_name} does not have enough columns.")
            continue

        # Extract sample name
        sample = '.'.join(file_name.split('.')[:2])

        # Calculate genomic inflation (lambda)
        x = chi2.ppf(1 - df[3], 1)
        lambda_val = np.median(x) / chi2.ppf(0.5, 1)

        # Check if lambda is within acceptable range
        if not (0.975 <= lambda_val <= 1.025):
            remove_by_lambda.append(sample)

    except Exception as e:
        print(f"Error processing file {file_name}: {e}")
        continue

# Write results to the output file
output_file = os.path.join(out_dir, "01.0_nonSig_by_lambda")
with open(output_file, "w") as f:
    for sample in remove_by_lambda:
        f.write(sample + '\n')

print(f"Samples failing lambda check written to: {output_file}")

