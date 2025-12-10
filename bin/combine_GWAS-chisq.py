import os
import gzip
import numpy as np
from scipy.stats import combine_pvalues
import argparse

def process_chunk(chunk, input_files, key_index_map, chunk_size):
    pvalue_array = np.full((chunk_size, len(input_files)), 0.9999)
    file_counts = np.zeros(chunk_size, dtype=int)
    keys = [None] * chunk_size

    for file_idx, file in enumerate(input_files):
        with gzip.open(file, 'rt') as f:
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
                            file_counts[row_idx] += 1
                            keys[row_idx] = key
                except (ValueError, IndexError):
                    continue

    return keys, pvalue_array, file_counts

def combine_chisq_for_chunk(pvalue_array):
    combined_chisq = np.apply_along_axis(lambda row: combine_pvalues(row, method='fisher')[0], 1, pvalue_array)
    return combined_chisq

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

    with open(os.path.join(output_dir, f'combined_chisq_{dmr_type}.csv'), 'w') as outfile:
        outfile.write("key\tcombined_chisq\tfile_count\n")
        for chunk in range(num_chunks):
            print(f"Processing chunk {chunk + 1}/{num_chunks}...")
            keys, pvalue_array, file_counts = process_chunk(chunk, input_files, key_index_map, chunk_size)
            combined_chisq = combine_chisq_for_chunk(pvalue_array)
            for i in range(chunk_size):
                if keys[i] is not None:
                    outfile.write(f"{keys[i]}\t{combined_chisq[i]:.6f}\t{file_counts[i]}\n")

    print("Processing complete!")

if __name__ == '__main__':
    main()
