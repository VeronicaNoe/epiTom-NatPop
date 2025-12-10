import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Define chromosome sizes
chr_sizes = {
    '01': 98543444,
    '02': 55340444,
    '03': 70787664,
    '04': 66470942,
    '05': 65875088,
    '06': 49751636,
    '07': 68045021,
    '08': 65866657,
    '09': 72482091,
    '10': 65527505,
    '11': 56302525,
    '12': 67145203,
}

def process_data(file_prefix):
    # Load the data
    df = pd.read_csv(f"02.1_combined_pvalues_{file_prefix}.csv", sep="\t")
    df['combined_pvalue'] = pd.to_numeric(df['combined_pvalue'], errors='coerce')
    # Split the 'key' column into chromosome and position
    df[['CHR', 'BP']] = df['key'].str.split(':', expand=True)
    df['BP'] = df['BP'].astype(int)
    df['CHR'] = df['CHR'].astype(str).str.zfill(2)
    # Create cumulative base pair column
    current_position = 0
    chr_positions = {}
    for chr, size in chr_sizes.items():
        chr_positions[chr] = current_position
        current_position += size + 1  # Add 1 bp to start the next chromosome
    df['cumulative_BP'] = df.apply(lambda row: row['BP'] + chr_positions[row['CHR']], axis=1)
    # Calculate top 1% combined_chisq
    top_1_percent = int(len(df) * 0.01)
    df_sorted = df.sort_values(by='combined_chisq', ascending=False)  # Corrected line
    df_top_1_percent = df_sorted.head(top_1_percent)
    df_top_1_percent = df_top_1_percent.sort_values(by=['CHR', 'BP'], ascending=[True, True])
    df_top_1_percent.to_csv(f"02.2_{file_prefix}_top_1_percent.csv", index=False)
    df_top_1_percent['key'].to_csv(f"02.2_{file_prefix}_top_1_percent_keys.txt", index=False, header=False)

# Process both C-DMR and CG-DMR files
for prefix in ['C-DMR', 'CG-DMR']:
    process_data(prefix)
