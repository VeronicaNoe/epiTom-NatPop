import os
import numpy as np
from scipy.stats import combine_pvalues
import gzip
import argparse
from multiprocessing import cpu_count, Pool, Manager
import warnings

def process_line(line):
    parts = line.strip().split()
    try:
        key = parts[0]  # Use the first column as the key
        pvalue = float(parts[3])  # Use the fourth column as the p-value
        return key, pvalue
    except (ValueError, IndexError):
        return None


def read_file_chunk(file, chunk_size, start_line):
    combined = {}
    with gzip.open(file, 'rt') as f:
        for line_number in range(start_line):
            next(f, None)  # Skip lines to reach the start of the chunk
        for line_number in range(chunk_size):
            line = f.readline()
            if not line:
                break
            result = process_line(line)
            if result is not None:
                key, pvalue = result
                if key not in combined:
                    combined[key] = []
                combined[key].append(pvalue)
    return combined
    


def combine_pvalues_fisher(pvalues):
    pvalues = np.where(np.isnan(pvalues), 0.999, pvalues)
    if not pvalues:
        return np.nan, np.nan, 0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        combined_chisq = -2 * np.sum(np.log(pvalues))
        combined_pvalue = combine_pvalues(pvalues, method='fisher')[1]
    
    return combined_pvalue, combined_chisq, len(pvalues)


def process_chunk(chunk_num, input_files, chunk_size):
    start_line = chunk_size * chunk_num
    chunk_results = {}
    
    for file in input_files:
        combined = read_file_chunk(file, chunk_size, start_line)
        for key, pvalues in combined.items():
            if key not in chunk_results:
                chunk_results[key] = pvalues
            else:
                chunk_results[key].extend(pvalues)

    combined_chunk_results = {}
    for key, pvalues in chunk_results.items():
        combined_pvalue, combined_chisq, pvalues_len = combine_pvalues_fisher(pvalues)
        combined_chunk_results[key] = [combined_chisq, combined_pvalue, pvalues_len]
    
    return combined_chunk_results

def accumulate_results(input_files, chunk_size, num_cpus):
    manager = Manager()
    global_results = manager.dict()
    max_lines = 1300000  # Maximum number of lines in any file
    num_chunks = (max_lines // chunk_size) + 1

    with Pool(processes=num_cpus) as pool:
        chunk_results = pool.starmap(process_chunk, [(chunk_num, input_files, chunk_size) for chunk_num in range(num_chunks)])
        
    for chunk_result in chunk_results:
        for key, values in chunk_result.items():
            if key not in global_results:
                global_results[key] = values
            else:
                global_results[key][0] += values[0]  # Sum combined_chisq
                global_results[key][1] = combine_pvalues([global_results[key][1], values[1]], method='fisher')[1]  # Combine p-values
                global_results[key][2] += values[2]  # Sum pvalues_len
    
    return dict(global_results)

def merge_results(global_results):
    final_results = []
    for key, values in global_results.items():
        combined_chisq, combined_pvalue, pvalues_len = values
        final_results.append([key, combined_chisq, combined_pvalue, min(global_results[key]), pvalues_len])
    return final_results

def main():
    parser = argparse.ArgumentParser(description='Process DMR files.')
    parser.add_argument('dmr_type', type=str, help='Type of DMR files to process (e.g., C-DMR, CG-DMR)')
    parser.add_argument('--num_files', type=int, default=None, help='Number of files to process (default: all)')
    parser.add_argument('--chunk_size', type=int, default=100000, help='Number of lines to read per chunk (default: 100000)')
    parser.add_argument('--num_cpus', type=int, default=cpu_count(), help='Number of CPUs to use (default: all available CPUs)')
    args = parser.parse_args()

    dmr_type = args.dmr_type
    num_files = args.num_files
    chunk_size = args.chunk_size
    num_cpus = args.num_cpus

    input_dir = "/mnt/disk3/vibanez/DMR-GWAS/data"  # Update to your data directory
    output_dir = "/mnt/disk3/vibanez/DMR-GWAS/results"
    
    os.chdir(input_dir)
    input_files = [file for file in os.listdir() if file.endswith('.ps.gz') and dmr_type in file]

    if num_files is not None:
        input_files = input_files[:num_files]

    global_results = accumulate_results(input_files, chunk_size, num_cpus)

    final_combined_results = merge_results(global_results)

    output_file = os.path.join(output_dir, f'combined_results_{dmr_type}.tsv')
    with open(output_file, 'w') as f:
        f.write('V1\tcombined_chisq\tcombined_pvalue\tmin_p_value\tnDMR\n')
        for row in final_combined_results:
            f.write('\t'.join(map(str, row)) + '\n')

if __name__ == '__main__':
    main()
