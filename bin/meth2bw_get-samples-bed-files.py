import sys
import os

INDIR = "/mnt/disk2/vibanez/02_methylkit/ac_annotation/ca_mergedAnnotation"
OUTDIR = "/mnt/disk2/vibanez/02_methylkit/ac_annotation/be_bigwig-files/bed-files"

def split_file(main_file_name, sample_names_file):
    # Create the full path for the main file
    main_file = os.path.join(INDIR, main_file_name)
    # Read sample names from the sample file
    with open(sample_names_file, 'r') as sample_file:
        sample_names = [line.strip() for line in sample_file]
    # Create a dictionary to store file handles for each sample
        sample_files = {sample_name: open(os.path.join(OUTDIR, f"{sample_name}.{main_file_name}.tmp"), 'w') for sample_name in sample_names}
    # Read the main file and split data into separate files based on sample names
    with open(main_file, 'r') as main_file:
        for line in main_file:
            # Split the line into columns
            columns = line.strip().split('\t')
            region_info = '\t'.join(columns[:3])
            # Write the region info to each sample's file
            for sample_name in sample_names:
                sample_data = '\t'.join([region_info] + [columns[3 + sample_names.index(sample_name)]])
                sample_files[sample_name].write(sample_data + '\n')
    # Close all output files
    for file_handle in sample_files.values():
        file_handle.close()


if __name__ == "__main__":
    main_file_name = sys.argv[1]  # Path to the main file
    sample_names_file = "00_colNames.tsv"  # Path to the sample names file
    split_file(main_file_name, sample_names_file)
