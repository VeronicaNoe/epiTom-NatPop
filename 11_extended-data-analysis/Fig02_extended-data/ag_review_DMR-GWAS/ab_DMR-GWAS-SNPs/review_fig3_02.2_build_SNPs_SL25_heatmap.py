#!/users/bioinfo/vibanez/anaconda3/bin/python

import os
import re
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, save_npz

# =========================
# Arguments
# =========================
parser = argparse.ArgumentParser()
parser.add_argument("dmr_type", choices=["C-DMR", "CG-DMR"])
parser.add_argument("--bin-size", type=int, default=1_000_000)
args = parser.parse_args()

DMR_TYPE = args.dmr_type
bin_size = args.bin_size
P_FLOOR = 1e-300
PRIMARY_CHROMS = [f"{i:02d}" for i in range(1, 13)]

# =========================
# Paths
# =========================
BASE_INPUT_DIR = "/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results/sig"
RESULTS_DIR = "/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/ab_DMR-GWAS-SNPs"
PLOTS_DIR = "/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/ab_DMR-GWAS-SNPs"

BIM_FILE = "/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/ba_markers/SNPs.LD.bim"
SL25_FAI = "/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/ITAG2.4_genomic.fasta.fai"

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

FINAL_MATRIX_FILE = os.path.join(RESULTS_DIR, f"heatmap_{DMR_TYPE}_SNP_SL25_binned.npz")

# =========================
# Helpers
# =========================
def normalize_chr(chrn):
    s = str(chrn).strip()
    s = re.sub(r"^SL2\.50ch", "", s, flags=re.IGNORECASE)
    s = re.sub(r"^(chr|ch)", "", s, flags=re.IGNORECASE)

    m = re.search(r"(\d+)$", s)
    if m is None:
        return None

    return f"{int(m.group(1)):02d}"


def load_chr_sizes_from_fai(fai_path):
    chr_sizes = {}
    with open(fai_path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            chrom = normalize_chr(parts[0])
            if chrom is None:
                continue

            if chrom in PRIMARY_CHROMS:
                chr_sizes[chrom] = int(parts[1])

    missing = [c for c in PRIMARY_CHROMS if c not in chr_sizes]
    if missing:
        raise ValueError(f"Missing chromosomes in FAI: {missing}")

    return chr_sizes


def linear_genome_pos(chrn, pos):
    return cumulative_lengths[chrn] + pos


def dmr_id_from_filename(file_path):
    base = os.path.basename(file_path)
    return re.sub(r"\.mQTL$", "", base)


def parse_dmr_sl25_from_id(dmr_id):
    parts = dmr_id.split("_")
    if len(parts) != 3:
        return None

    chrn = normalize_chr(parts[0])
    if chrn is None or chrn not in PRIMARY_CHROMS:
        return None

    try:
        start = int(parts[2])
    except ValueError:
        return None

    end = start + 99
    mid = (start + end) // 2
    return chrn, start, end, mid


def load_snp_coords_from_bim(path, chr_sizes):
    coord_map = {}
    with open(path, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue

            snp_id = parts[1]
            chrom = normalize_chr(parts[0])
            if chrom is None or chrom not in PRIMARY_CHROMS:
                continue

            try:
                pos = int(parts[3])
            except ValueError:
                continue

            if pos < 0 or pos > chr_sizes[chrom]:
                continue

            coord_map[snp_id] = (chrom, pos)

    return coord_map


# =========================
# Load genome layout
# =========================
chr_sizes = load_chr_sizes_from_fai(SL25_FAI)
ordered_chroms = sorted(chr_sizes.keys(), key=int)

print(f"Loaded SL2.5 chromosome sizes from: {SL25_FAI}")
for chrom in ordered_chroms:
    print(f"{chrom}\t{chr_sizes[chrom]}")

cumulative_lengths = {}
cum = 0
for chrom in ordered_chroms:
    cumulative_lengths[chrom] = cum
    cum += chr_sizes[chrom]

genome_length = cum
num_bins = genome_length // bin_size + 1

print(f"Genome length: {genome_length}")
print(f"Number of bins: {num_bins}")

heatmap_data_binned = lil_matrix((num_bins, num_bins), dtype=np.float32)

# =========================
# Load SNP coordinates
# =========================
snp_coord_map = load_snp_coords_from_bim(BIM_FILE, chr_sizes)
print(f"Loaded SNP coordinates: {len(snp_coord_map)}")

# =========================
# Collect input files
# =========================
input_files = sorted(glob.glob(os.path.join(BASE_INPUT_DIR, f"*{DMR_TYPE}*.mQTL")))
print(f"Found input files: {len(input_files)}")

# =========================
# Build matrix
# =========================
missing_dmr = 0
missing_snp = 0
kept_lines = 0

for i, file_path in enumerate(input_files, start=1):
    dmr_id = dmr_id_from_filename(file_path)
    parsed = parse_dmr_sl25_from_id(dmr_id)
    if parsed is None:
        missing_dmr += 1
        continue

    y_chr, y_start, y_end, y_mid = parsed
    y_idx = linear_genome_pos(y_chr, y_mid) // bin_size

    with open(file_path, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue

            snp_id = parts[0]
            try:
                pval = float(parts[3])
            except ValueError:
                continue

            if np.isnan(pval):
                continue
            if pval <= 0:
                pval = P_FLOOR

            if snp_id not in snp_coord_map:
                missing_snp += 1
                continue

            chrn, pos = snp_coord_map[snp_id]
            x_idx = linear_genome_pos(chrn, pos) // bin_size

            current = heatmap_data_binned[y_idx, x_idx]
            if current == 0 or pval < current:
                heatmap_data_binned[y_idx, x_idx] = pval

            kept_lines += 1

    if i % 100 == 0:
        print(f"Processed {i}/{len(input_files)} files")

heatmap_data_binned = heatmap_data_binned.tocsr()
save_npz(FINAL_MATRIX_FILE, heatmap_data_binned)

print(f"Saved: {FINAL_MATRIX_FILE}")
print(f"Missing DMR IDs: {missing_dmr}")
print(f"Missing SNP IDs: {missing_snp}")
print(f"Stored signals: {heatmap_data_binned.nnz}")
