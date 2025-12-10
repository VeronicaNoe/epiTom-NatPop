import os
import gzip
import numpy as np
from scipy.stats import combine_pvalues
import argparse
from multiprocessing import Pool, cpu_count


def process_chunk(chunk, input_files, key_index_map, chunk_size):
    # Initialize the array for this chunk with 0.9999
    pvalue_array = np.full((chunk_size, len(input_files)), 0.9999)
    file_counts = np.zeros(chunk_size, dtype=int)

    keys = [None] * chunk_size  # To store the key for each row

    for file_idx, file in enumerate(input_files):
        with gzip.open(file, 'rt') as f:
            # Skip lines to reach the start of the chunk
            for _ in range(chunk * chunk_size):
                next(f, None)
            for line_idx in range(chunk_size):
                line = f.readline()
                if not line:
                    break
                parts = line.strip().split()
                try:
                    key = parts[0]
                    pvalue = float(parts[3]) if parts[3] != 'nan' else 0.9999
                    if key in key_index_map:
                        row_idx = key_index_map[key] - (chunk * chunk_size)
                        if 0 <= row_idx < chunk_size:
                            pvalue_array[row_idx, file_idx] = pvalue
                            file_counts[row_idx] += 1  # Count the number of files contributing to each row
                            keys[row_idx] = key  # Store the key for this row
                except (ValueError, IndexError):
                    continue

    return keys, pvalue_array, file_counts

def combine_pvalues_for_chunk(pvalue_array):
    combined_results = np.apply_along_axis(lambda row: combine_pvalues(row, method='fisher'), 1, pvalue_array)
    combined_chisq = combined_results[:, 0]
    combined_pvalues = combined_results[:, 1]
    return combined_chisq, combined_pvalues

def main():
    parser = argparse.ArgumentParser(description='Process DMR files.')
    parser.add_argument('dmr_type', type=str, help='Type of DMR files to process (e.g., C-DMR, CG-DMR)')
    args = parser.parse_args()

    dmr_type = args.dmr_type
    input_dir = "/mnt/disk3/vibanez/DMR-GWAS/data"
    output_dir = "/mnt/disk3/vibanez/DMR-GWAS/results"

    os.chdir(input_dir)
    input_files = [file for file in os.listdir() if file.endswith('.ps.gz') and dmr_type in file]

    if not input_files:
        print("Error: No input files found.")
        return

    total_keys = 1228690
    chunk_size = 10000

    key_index_map = {}
    with gzip.open(input_files[0], 'rt') as file:
        for idx, line in enumerate(file):
            key = line.split()[0]
            key_index_map[key] = idx

    num_chunks = (total_keys // chunk_size) + 1

    # Use all available CPU cores
    num_workers = cpu_count()
    
    # Create a multiprocessing pool
    with Pool(num_workers) as pool:
        # Parallelize the processing of chunks
        results = pool.starmap(process_chunk, [(chunk, input_files, key_index_map, chunk_size) for chunk in range(num_chunks)])

    with open(os.path.join(output_dir, f'combined_pvalues_{dmr_type}.csv'), 'w') as outfile:
        outfile.write("key\tcombined_chisq\tcombined_pvalue\tfile_count\n")
        for chunk_idx, (keys, pvalue_array, file_counts) in enumerate(results):
            print(f"Processing chunk {chunk_idx + 1}/{num_chunks}...")
            combined_chisq, combined_pvalues = combine_pvalues_for_chunk(pvalue_array)
            for i in range(chunk_size):
                if keys[i] is not None:
                    outfile.write(f"{keys[i]}\t{combined_chisq[i]:.6f}\t{combined_pvalues[i]:.6e}\t{file_counts[i]}\n")

    print("Processing complete!")

if __name__ == '__main__':
    main()
