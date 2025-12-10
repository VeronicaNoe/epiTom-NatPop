#!/usr/bin/env python3

import argparse
import pandas as pd

def main(ld_file, output_file, max_dist, bin_size):
    # Read PLINK LD file
    df = pd.read_csv(ld_file, sep="\s+", comment="#")

    # Calculate midpoint and mid-distance
    df["mid_dist"] = abs((df["BP_B"] + df["BP_A"]) / 2 - df["BP_A"])
    df["midpoint"] = (df["BP_A"] + df["BP_B"]) // 2

    # Filter by max distance
    df = df[df["mid_dist"] <= max_dist]

    # Assign bins
    df["bin"] = (df["mid_dist"] // bin_size).astype(int) * bin_size

    # Summarize
    summary = (
        df.groupby("bin")
        .agg(
            mean_r2=("R2", "mean"),
            median_r2=("R2", "median"),
            n_pairs=("R2", "count"),
            mean_midpoint=("midpoint", "mean")
        )
        .reset_index()
        .sort_values("bin")
    )

    # Add bin edges
    summary["bin_start"] = summary["bin"]
    summary["bin_end"] = summary["bin_start"] + bin_size
    summary = summary[["bin_start", "bin_end", "mean_r2", "median_r2", "n_pairs", "mean_midpoint"]]

    # Save
    summary.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Estimate LD decay from PLINK r2 file')
    parser.add_argument('--ld', required=True, help='plink ld.gz file (r2 output)')
    parser.add_argument('--out', required=True, help='output file')
    parser.add_argument('--max_dist', type=int, default=100000, help='max distance to consider')
    parser.add_argument('--step', type=int, default=100, help='bin size for distance')
    args = parser.parse_args()

    main(ld_file=args.ld, output_file=args.out, max_dist=args.max_dist, bin_size=args.step)
