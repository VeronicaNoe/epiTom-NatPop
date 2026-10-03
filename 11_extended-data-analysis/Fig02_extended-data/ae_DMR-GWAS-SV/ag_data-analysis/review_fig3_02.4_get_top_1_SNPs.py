import pandas as pd
import numpy as np

def process_data(file_prefix):
    in_file = f"/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/results/02.1_combined_pvalues_{file_prefix}.csv"
    out_csv = f"/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/results/02.2_{file_prefix}_top_1_percent.csv"
    out_keys = f"/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/results/02.2_{file_prefix}_top_1_percent_keys.txt"

    # Load data
    df = pd.read_csv(in_file, sep="\t")

    # Use coordinates already present in the file
    df["CHR"] = df["CHR"].astype(str).str.zfill(2)
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    df["combined_pvalue"] = pd.to_numeric(df["combined_pvalue"], errors="coerce")
    df["combined_chisq"] = pd.to_numeric(df["combined_chisq"], errors="coerce")

    # Keep valid rows
    df = df.dropna(subset=["CHR", "BP", "combined_chisq"]).copy()
    df["BP"] = df["BP"].astype(int)

    # Top 1%
    top_1_percent = max(1, int(len(df) * 0.01))
    df_sorted = df.sort_values(by="combined_chisq", ascending=False)
    df_top_1_percent = df_sorted.head(top_1_percent).copy()

    # Sort output by genomic coordinate
    df_top_1_percent = df_top_1_percent.sort_values(by=["CHR", "BP"], ascending=[True, True])

    # Save full table with coordinates
    df_top_1_percent.to_csv(out_csv, index=False)

    # Save only keys
    df_top_1_percent["key"].to_csv(out_keys, index=False, header=False)

for prefix in ["C-DMR", "CG-DMR"]:
    process_data(prefix)
