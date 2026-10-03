#!/users/bioinfo/vibanez/anaconda3/bin/python

import os
import re
import gzip
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.sparse import lil_matrix, save_npz


# =========================
# Arguments
# =========================
parser = argparse.ArgumentParser()
parser.add_argument("dmr_type", choices=["C-DMR", "CG-DMR"], help="DMR type to process")
parser.add_argument("--threshold", type=float, default=6.492663e-06, help="P-value threshold for heatmap")
parser.add_argument("--bin-size", type=int, default=1_000_000, help="Bin size in bp")
args = parser.parse_args()

DMR_TYPE = args.dmr_type
threshold = args.threshold
bin_size = args.bin_size
P_FLOOR = 1e-300

PRIMARY_CHROMS = [f"{i:02d}" for i in range(1, 13)]


# =========================
# Paths
# =========================
BASE_INPUT_DIR = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ad_results"
BASE_ANALYSIS_DIR = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis"
RESULTS_DIR = "/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS"
PLOTS_DIR = "/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS"

SV_COORD_FILE = "/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/sv_id2coord_sl25_primary.tsv"
SL25_FAI = "/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/ITAG2.4_genomic.fasta.fai"

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

COMBINED_TSV = os.path.join(RESULTS_DIR, f"combined_fisher_SV_{DMR_TYPE}.tsv")
FINAL_MATRIX_FILE = os.path.join(RESULTS_DIR, f"heatmap_{DMR_TYPE}_SV_SL25_binned.npz")


# =========================
# Helpers
# =========================
def normalize_chr(chrn):
    s = str(chrn).strip()
    s = re.sub(r"^(SL2\.50ch|chr|ch)", "", s, flags=re.IGNORECASE)

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


def linear_genome_pos(chrn, pos):
    return cumulative_lengths[chrn] + pos


def dmr_id_from_filename(file_path):
    base = os.path.basename(file_path)
    return re.sub(r"\.ps\.gz$", "", base)


def parse_dmr_sl25_from_id(dmr_id):
    # example: ch01_C-DMR_10001
    parts = dmr_id.split("_")
    if len(parts) != 3:
        return None

    chr_raw = parts[0]
    pos_raw = parts[2]

    chrn = normalize_chr(chr_raw)
    if chrn is None or chrn not in PRIMARY_CHROMS:
        return None

    try:
        start = int(pos_raw)
    except ValueError:
        return None

    end = start + 99
    mid = (start + end) // 2

    return chrn, start, end, mid


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


def save_combined_results(tsv_file, sv_ids, sum_log_p, file_count, sv_coord_map):
    valid_idx = np.where(file_count > 0)[0]

    with open(tsv_file, "w") as out:
        out.write("key\tchrom\tpos\tlinear_pos\tx_bin\tcombined_chisq\tcombined_pvalue\tfile_count\n")

        for i in valid_idx:
            sv_id = sv_ids[i]
            if sv_id not in sv_coord_map:
                continue

            chrom, pos = sv_coord_map[sv_id]
            linpos = linear_genome_pos(chrom, pos)
            x_bin = linpos / bin_size

            chisq = -2.0 * sum_log_p[i]
            pval = chi2.sf(chisq, 2 * file_count[i])

            out.write(
                f"{sv_id}\t{chrom}\t{pos}\t{linpos}\t{x_bin:.6f}\t"
                f"{chisq:.6f}\t{pval:.6e}\t{file_count[i]}\n"
            )

    print(f"Saved combined Fisher table: {tsv_file}")


# =========================
# Load chromosome sizes from SL2.5 FAI
# =========================
chr_sizes = load_chr_sizes_from_fai(SL25_FAI)
ordered_chroms = sorted(chr_sizes.keys(), key=int)

print(f"Loaded SL2.5 chromosome sizes from: {SL25_FAI}")
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

print(f"Genome length: {genome_length}")
print(f"Number of bins: {num_bins}")

# rows = DMR bins, cols = SV bins
heatmap_data_binned = lil_matrix((num_bins, num_bins), dtype=np.float32)

# =========================
# Load SV coordinates
# =========================
sv_coord_map = load_sv_coords(SV_COORD_FILE, chr_sizes)

# create stable SV order for combined Fisher accumulation
sv_ids = sorted(sv_coord_map.keys())
sv_index = {sv: i for i, sv in enumerate(sv_ids)}
sum_log_p = np.zeros(len(sv_ids), dtype=np.float64)
file_count = np.zeros(len(sv_ids), dtype=np.int32)

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
# Build heatmap + Fisher stats
# =========================
total_missing_sv = 0
total_kept_lines = 0
missing_dmr = 0
dmr_files_per_ychrom = {c: 0 for c in ordered_chroms}

for i, file_path in enumerate(input_files, start=1):
    dmr_id = dmr_id_from_filename(file_path)
    parsed = parse_dmr_sl25_from_id(dmr_id)

    if parsed is None:
        missing_dmr += 1
        continue

    y_chr, y_start, y_end, y_mid = parsed
    y_linear = linear_genome_pos(y_chr, y_mid)
    y_idx = y_linear // bin_size

    if y_idx < 0 or y_idx >= num_bins:
        missing_dmr += 1
        continue

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

            # ---- combined Fisher: all valid p-values
            idx = sv_index.get(sv_id)
            if idx is not None:
                sum_log_p[idx] += np.log(pval)
                file_count[idx] += 1

            # ---- heatmap: only significant values
            if pval >= threshold:
                continue

            if sv_id not in sv_coord_map:
                total_missing_sv += 1
                continue

            chrn, pos = sv_coord_map[sv_id]
            linear_pos = linear_genome_pos(chrn, pos)
            x_idx = linear_pos // bin_size

            if x_idx < 0 or x_idx >= num_bins:
                total_missing_sv += 1
                continue

            current = heatmap_data_binned[y_idx, x_idx]
            if current == 0 or pval < current:
                heatmap_data_binned[y_idx, x_idx] = pval

            total_kept_lines += 1

    dmr_files_per_ychrom[y_chr] += 1

    if i % 100 == 0:
        print(f"Processed {i}/{len(input_files)} files")
        print(f"  missing DMRs so far: {missing_dmr}")
        print(f"  missing SV IDs so far: {total_missing_sv}")
        print(f"  kept significant lines so far: {total_kept_lines}")

        save_npz(
            os.path.join(RESULTS_DIR, f"partial_{DMR_TYPE}_SV_SL25_binned.npz"),
            heatmap_data_binned.tocsr()
        )

# =========================
# Save outputs
# =========================
heatmap_data_binned = heatmap_data_binned.tocsr()
save_npz(FINAL_MATRIX_FILE, heatmap_data_binned)

save_combined_results(COMBINED_TSV, sv_ids, sum_log_p, file_count, sv_coord_map)

print(f"Saved matrix: {FINAL_MATRIX_FILE}")
print(f"Total missing DMR IDs: {missing_dmr}")
print(f"Total missing SV IDs: {total_missing_sv}")
print(f"Total kept significant lines: {total_kept_lines}")
print(f"Total stored signals (nnz): {heatmap_data_binned.nnz}")

print("DMR files contributing per SL2.5 chromosome:")
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
    plt.xlabel("SV genome bin in SL2.5 (chr01-chr12, 1 Mb)")
    plt.ylabel(f"DMR genome bin in SL2.5 (chr{chrom}, 1 Mb)")
    plt.title(f"{DMR_TYPE} SV-vs-DMR heatmap in SL2.5")

    output_file = os.path.join(PLOTS_DIR, f"{DMR_TYPE}_SV_SL25_genome_heatmap_chr{chrom}.pdf")
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Saved: {output_file}")
