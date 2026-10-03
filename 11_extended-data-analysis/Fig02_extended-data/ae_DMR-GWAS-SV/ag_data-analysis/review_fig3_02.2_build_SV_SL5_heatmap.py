#!/users/bioinfo/vibanez/anaconda3/bin/python

import os
import re
import gzip
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, save_npz, load_npz


# =========================
# Arguments
# =========================
parser = argparse.ArgumentParser()
parser.add_argument("dmr_type", choices=["C-DMR", "CG-DMR"], help="DMR type to process")
parser.add_argument("--threshold", type=float, default=6.492663e-06, help="P-value threshold")
parser.add_argument("--bin-size", type=int, default=1_000_000, help="Bin size in bp")
args = parser.parse_args()

DMR_TYPE = args.dmr_type
threshold = args.threshold
bin_size = args.bin_size
P_FLOOR = 1e-300


# =========================
# Keep only chr01-chr12
# =========================
PRIMARY_CHROMS = [f"{i:02d}" for i in range(1, 13)]


# =========================
# Paths
# =========================
BASE_INPUT_DIR = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ad_results"
BASE_ANALYSIS_DIR = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis"
RESULTS_DIR = os.path.join(BASE_ANALYSIS_DIR, "results")
PLOTS_DIR = os.path.join(BASE_ANALYSIS_DIR, "plots")

LIFTOFF_DIR = "/mnt/disk2/vibanez/otherAnalysis/01_liftoff-DMRs-SL2.5-SL5"
DMR_COORD_FILE = os.path.join(LIFTOFF_DIR, f"{DMR_TYPE}_id_to_SL5.tsv")

SV_COORD_FILE = os.path.join(BASE_ANALYSIS_DIR, "bb_merge-results", "sv.id2coord.tsv")
SL5_FAI = "/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/SL5.0.fasta.fai"

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)


# =========================
# Helpers
# =========================
def normalize_chr(chrn):
    s = str(chrn).strip()
    s = re.sub(r"^(chr|ch)", "", s, flags=re.IGNORECASE)

    m = re.search(r"(\d+)$", s)
    if m is None:
        return None

    return f"{int(m.group(1)):02d}"


def linear_genome_pos(chrn, pos):
    return cumulative_lengths[chrn] + pos


def dmr_id_from_filename(file_path):
    base = os.path.basename(file_path)
    return re.sub(r"\.ps\.gz$", "", base)


def load_chr_sizes_from_fai(fai_path):
    chr_sizes = {}

    with open(fai_path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            chrom_raw = parts[0]
            size = int(parts[1])

            chrom = normalize_chr(chrom_raw)
            if chrom is None:
                continue

            if chrom in PRIMARY_CHROMS:
                chr_sizes[chrom] = size

    missing = [c for c in PRIMARY_CHROMS if c not in chr_sizes]
    if missing:
        raise ValueError(f"Missing primary chromosomes in FAI: {missing}")

    return chr_sizes


def load_dmr_sl5_coords_from_tsv(tsv_path, dmr_type, chr_sizes):
    """
    Expected TSV format like:
    SL2.50ch01:901-1000    1    85378    85478

    We convert that to:
    ch01_C-DMR_901  ->  (sl5_chr, sl5_start, sl5_end, sl5_mid)
    """
    dmr_map = {}
    duplicate_count = 0
    skipped_bad_chr = 0
    skipped_bad_pos = 0
    skipped_bad_id = 0

    with open(tsv_path, "rt") as f:
        for line in f:
            if not line.strip():
                continue

            parts = line.strip().split()
            if len(parts) < 4:
                continue

            src_id = parts[0]
            sl5_chr_raw = parts[1]
            sl5_start_raw = parts[2]
            sl5_end_raw = parts[3]

            m = re.match(r"SL2\.50ch(\d+):(\d+)-(\d+)$", src_id)
            if m is None:
                skipped_bad_id += 1
                continue

            src_chr = m.group(1).zfill(2)
            src_start = int(m.group(2))

            dmr_id = f"ch{src_chr}_{dmr_type}_{src_start}"

            sl5_chr = normalize_chr(sl5_chr_raw)
            if sl5_chr is None or sl5_chr not in PRIMARY_CHROMS:
                skipped_bad_chr += 1
                continue

            if sl5_chr not in chr_sizes:
                skipped_bad_chr += 1
                continue

            try:
                sl5_start = int(sl5_start_raw)
                sl5_end = int(sl5_end_raw)
            except ValueError:
                skipped_bad_pos += 1
                continue

            if sl5_start > sl5_end:
                sl5_start, sl5_end = sl5_end, sl5_start

            if sl5_start < 0 or sl5_end > chr_sizes[sl5_chr]:
                skipped_bad_pos += 1
                continue

            sl5_mid = (sl5_start + sl5_end) // 2

            if dmr_id in dmr_map:
                duplicate_count += 1
                continue

            dmr_map[dmr_id] = (sl5_chr, sl5_start, sl5_end, sl5_mid)

    print(f"Loaded {len(dmr_map)} lifted DMR coordinates from: {tsv_path}")
    print(f"Duplicate lifted DMR entries skipped: {duplicate_count}")
    print(f"Lifted DMRs skipped outside chr01-chr12: {skipped_bad_chr}")
    print(f"Lifted DMRs skipped due to bad position: {skipped_bad_pos}")
    print(f"Lifted DMRs skipped due to bad source ID: {skipped_bad_id}")

    return dmr_map


def load_sv_coords(path, chr_sizes):
    sv_coord_map = {}
    skipped_bad_chr = 0
    skipped_bad_pos = 0

    with open(path, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue

            sv_id = parts[0]
            chrom_raw = parts[1]
            pos_raw = parts[2]

            chrom = normalize_chr(chrom_raw)
            if chrom is None or chrom not in PRIMARY_CHROMS:
                skipped_bad_chr += 1
                continue

            if chrom not in chr_sizes:
                skipped_bad_chr += 1
                continue

            try:
                pos = int(pos_raw)
            except ValueError:
                skipped_bad_pos += 1
                continue

            if pos < 0 or pos > chr_sizes[chrom]:
                skipped_bad_pos += 1
                continue

            sv_coord_map[sv_id] = (chrom, pos)

    print(f"Loaded {len(sv_coord_map)} SV coordinates from: {path}")
    print(f"SVs skipped outside chr01-chr12: {skipped_bad_chr}")
    print(f"SVs skipped due to bad position: {skipped_bad_pos}")

    return sv_coord_map


def process_file(file_path, y_idx):
    missing_sv = 0
    kept_lines = 0

    with gzip.open(file_path, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue

            sv_id = parts[0]

            try:
                pval = float(parts[3])
            except ValueError:
                continue

            if np.isnan(pval):
                continue

            if pval <= 0:
                pval = P_FLOOR

            if pval >= threshold:
                continue

            if sv_id not in sv_coord_map:
                missing_sv += 1
                continue

            chrn, pos = sv_coord_map[sv_id]
            linear_pos = linear_genome_pos(chrn, pos)
            x_idx = linear_pos // bin_size

            if x_idx < 0 or x_idx >= num_bins:
                missing_sv += 1
                continue

            current = heatmap_data_binned[y_idx, x_idx]
            if current == 0 or pval < current:
                heatmap_data_binned[y_idx, x_idx] = pval

            kept_lines += 1

    return missing_sv, kept_lines


# =========================
# Load chromosome sizes from SL5 FAI
# =========================
chr_sizes = load_chr_sizes_from_fai(SL5_FAI)
ordered_chroms = sorted(chr_sizes.keys(), key=int)

print(f"Loaded SL5 chromosome sizes from: {SL5_FAI}")
for chrom in ordered_chroms:
    print(f"{chrom}\t{chr_sizes[chrom]}")


# =========================
# Linear genome coordinates
# =========================
cumulative_lengths = {}
cumulative_len = 0
for chr_id in ordered_chroms:
    cumulative_lengths[chr_id] = cumulative_len
    cumulative_len += chr_sizes[chr_id]

genome_length = cumulative_len
num_bins = genome_length // bin_size + 1

# rows = DMR bins, cols = SV bins
heatmap_data_binned = lil_matrix((num_bins, num_bins), dtype=np.float32)

print(f"Genome length: {genome_length}")
print(f"Number of bins: {num_bins}")


# =========================
# Load DMR and SV coordinates
# =========================
dmr_sl5_map = load_dmr_sl5_coords_from_tsv(DMR_COORD_FILE, DMR_TYPE, chr_sizes)
sv_coord_map = load_sv_coords(SV_COORD_FILE, chr_sizes)


# =========================
# Collect input files
# =========================
input_files = []
for chr_num in range(1, 13):
    chr_dir = f"ch{chr_num:02d}"
    pattern = os.path.join(BASE_INPUT_DIR, DMR_TYPE, chr_dir, "sig", "*.ps.gz")
    input_files.extend(sorted(glob.glob(pattern)))

print(f"Found {len(input_files)} files")


# =========================
# Build matrix
# =========================
total_missing_sv = 0
total_kept_lines = 0
missing_dmr_liftoff = 0

dmr_files_per_ychrom = {c: 0 for c in ordered_chroms}

for i, file_path in enumerate(input_files, start=1):
    dmr_id = dmr_id_from_filename(file_path)

    if dmr_id not in dmr_sl5_map:
        missing_dmr_liftoff += 1
        continue

    y_chr, y_start, y_end, y_mid = dmr_sl5_map[dmr_id]
    y_linear = linear_genome_pos(y_chr, y_mid)
    y_idx = y_linear // bin_size

    if y_idx < 0 or y_idx >= num_bins:
        missing_dmr_liftoff += 1
        continue

    missing_sv, kept_lines = process_file(file_path, y_idx)

    total_missing_sv += missing_sv
    total_kept_lines += kept_lines
    dmr_files_per_ychrom[y_chr] += 1

    if i % 100 == 0:
        print(f"Processed {i}/{len(input_files)} files")
        print(f"  missing lifted DMRs so far: {missing_dmr_liftoff}")
        print(f"  missing SV IDs so far: {total_missing_sv}")
        print(f"  kept significant lines so far: {total_kept_lines}")

        save_npz(
            os.path.join(RESULTS_DIR, f"partial_{DMR_TYPE}_SV_SL5_binned.npz"),
            heatmap_data_binned.tocsr()
        )


# =========================
# Save final matrix
# =========================
heatmap_data_binned = heatmap_data_binned.tocsr()
final_matrix_file = os.path.join(RESULTS_DIR, f"heatmap_{DMR_TYPE}_SV_SL5_binned.npz")
save_npz(final_matrix_file, heatmap_data_binned)

print(f"Saved matrix: {final_matrix_file}")
print(f"Total missing lifted DMRs: {missing_dmr_liftoff}")
print(f"Total missing SV IDs: {total_missing_sv}")
print(f"Total kept significant lines: {total_kept_lines}")
print(f"Total stored signals (nnz): {heatmap_data_binned.nnz}")

print("DMR files contributing per lifted SL5 chromosome:")
for chrom in ordered_chroms:
    print(f"{chrom}\t{dmr_files_per_ychrom[chrom]}")

print("Stored heatmap nnz per DMR chromosome block:")
for chrom in ordered_chroms:
    chr_start = cumulative_lengths[chrom]
    chr_end = chr_start + chr_sizes[chrom]
    chr_start_bin = chr_start // bin_size
    chr_end_bin = chr_end // bin_size + 1
    block_nnz = heatmap_data_binned[chr_start_bin:chr_end_bin, :].nnz
    print(f"{chrom}\t{block_nnz}")


# =========================
# Plot per chromosome
# =========================
for chrom in ordered_chroms:
    chr_start = cumulative_lengths[chrom]
    chr_end = chr_start + chr_sizes[chrom]

    chr_start_bin = chr_start // bin_size
    chr_end_bin = chr_end // bin_size + 1

    heatmap_data_chr = heatmap_data_binned[chr_start_bin:chr_end_bin, :]

    chr_rows, chr_cols = heatmap_data_chr.nonzero()
    chr_values = heatmap_data_chr.data

    if len(chr_values) == 0:
        print(f"No points for chromosome {chrom}")
        continue

    chr_values_transformed = -np.log10(chr_values)

    plt.figure(figsize=(15, 10))
    plt.scatter(
        chr_cols,
        chr_rows,
        c=chr_values_transformed,
        cmap="magma_r",
        marker="s",
        s=4,
        linewidths=0
    )

    for chrom_x in ordered_chroms:
        x_boundary = cumulative_lengths[chrom_x] // bin_size
        plt.axvline(x=x_boundary, color="grey", linestyle="--", linewidth=0.6, alpha=0.7)

    xticks = []
    xlabels = []
    for chrom_x in ordered_chroms:
        start = cumulative_lengths[chrom_x]
        end = start + chr_sizes[chrom_x]
        mid = ((start + end) // 2) // bin_size
        xticks.append(mid)
        xlabels.append(chrom_x)

    plt.xticks(xticks, xlabels, rotation=90)
    plt.colorbar(label="-log10(min p-value)")
    plt.xlabel("SV genome bin in SL5 (chr01-chr12, 1 Mb)")
    plt.ylabel(f"Lifted DMR genome bin in SL5 (chr{chrom}, 1 Mb)")
    plt.title(f"{DMR_TYPE} SV-vs-DMR heatmap in SL5")

    output_file = os.path.join(PLOTS_DIR, f"{DMR_TYPE}_SV_SL5_genome_heatmap_chr{chrom}.pdf")
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Saved: {output_file}")
