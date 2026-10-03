#!/users/bioinfo/vibanez/anaconda3/bin/python

import os
import re
import gzip
import argparse
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.stats import chi2
from scipy.sparse import load_npz


# =========================
# Arguments
# =========================
parser = argparse.ArgumentParser(
    description="Combine Fisher p-values across SV GWAS files and plot with corrected whole-genome heatmap."
)
parser.add_argument(
    "dmr_type",
    choices=["C-DMR", "CG-DMR"],
    help="Type of DMR files to process"
)
parser.add_argument(
    "--plot_only",
    action="store_true",
    help="Skip recomputing combined p-values and plot from existing TSV + NPZ"
)
args = parser.parse_args()

DMR_TYPE = args.dmr_type
PLOT_ONLY = args.plot_only


# =========================
# Paths
# =========================
BASE_DIR = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ad_results"
RESULTS_DIR = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/results"
PLOTS_DIR = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/plots"

SV_COORD_FILE = "/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ag_data-analysis/bb_merge-results/sv.id2coord.tsv"
SL5_FAI = "/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/SL5.0.fasta.fai"

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

COMBINED_TSV = os.path.join(RESULTS_DIR, f"combined_fisher_SV_{DMR_TYPE}.tsv")


# =========================
# Helpers
# =========================
def normalize_chr(chrn):
    chrn = str(chrn).strip()
    m = re.search(r"(\d+)$", chrn)
    if m is None:
        return None
    return m.group(1).zfill(2)


def load_chr_sizes_from_fai(fai_file):
    chr_sizes = {}

    with open(fai_file, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            raw_name = parts[0]
            size = int(parts[1])

            chrom = normalize_chr(raw_name)
            if chrom is None:
                continue

            if chrom in {f"{i:02d}" for i in range(1, 13)}:
                chr_sizes[chrom] = size

    chrom_order = [f"{i:02d}" for i in range(1, 13)]
    missing = [c for c in chrom_order if c not in chr_sizes]
    if missing:
        raise ValueError(f"Missing chromosomes in FAI: {missing}")

    return chr_sizes


def build_sv_pos_map(coord_file):
    sv_pos_map = {}

    with open(coord_file, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue

            sv_id = parts[0]
            chrom = normalize_chr(parts[1])
            if chrom is None:
                continue

            try:
                pos = int(parts[2])
            except ValueError:
                continue

            if chrom not in {f"{i:02d}" for i in range(1, 13)}:
                continue

            sv_pos_map[sv_id] = (chrom, pos)

    return sv_pos_map


def safe_pvalue(p):
    if not np.isfinite(p):
        return None
    if p <= 0:
        return 1e-300
    if p > 1:
        return None
    return p


def build_key_order(first_file):
    keys = []
    key_index_map = {}

    with gzip.open(first_file, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 1:
                continue

            key = parts[0]
            if key not in key_index_map:
                key_index_map[key] = len(keys)
                keys.append(key)

    return keys, key_index_map


def collect_input_files(base_dir, dmr_type, chrom_order):
    input_files = []
    for chrom in chrom_order:
        pattern = os.path.join(base_dir, dmr_type, f"ch{chrom}", "sig", "*.ps.gz")
        input_files.extend(sorted(glob(pattern)))
    return input_files


def linear_genome_pos(chrom, pos, cumulative_lengths):
    return cumulative_lengths[chrom] + pos


def save_combined_results(
    tsv_file,
    keys,
    combined_chisq,
    combined_pvalue,
    file_count,
    sv_pos_map,
    cumulative_lengths,
    bin_size
):
    n_written = 0

    with open(tsv_file, "w") as out:
        out.write("key\tchrom\tpos\tlinear_pos\tx_bin\tcombined_chisq\tcombined_pvalue\tfile_count\n")

        for i, key in enumerate(keys):
            if key not in sv_pos_map:
                continue
            if file_count[i] == 0:
                continue

            chrom, pos = sv_pos_map[key]
            if chrom not in cumulative_lengths:
                continue

            linpos = linear_genome_pos(chrom, pos, cumulative_lengths)
            x_bin = linpos / bin_size

            out.write(
                f"{key}\t{chrom}\t{pos}\t{linpos}\t{x_bin:.6f}\t"
                f"{combined_chisq[i]:.6f}\t{combined_pvalue[i]:.6e}\t{file_count[i]}\n"
            )
            n_written += 1

    print(f"Wrote {n_written} rows to: {tsv_file}")


def load_combined_results(tsv_file):
    chroms = []
    positions = []
    x_bins = []
    chisq = []
    pvals = []
    counts = []

    with open(tsv_file, "rt") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 8:
                continue

            try:
                chroms.append(parts[1])
                positions.append(int(parts[2]))
                x_bins.append(float(parts[4]))
                chisq.append(float(parts[5]))
                pvals.append(float(parts[6]))
                counts.append(int(parts[7]))
            except ValueError:
                continue

    return (
        np.array(chroms),
        np.array(positions, dtype=int),
        np.array(x_bins, dtype=float),
        np.array(chisq, dtype=float),
        np.array(pvals, dtype=float),
        np.array(counts, dtype=int),
    )


def find_heatmap_npz(results_dir, dmr_type):
    candidates = [
        os.path.join(results_dir, f"heatmap_{dmr_type}_SV_SL5_binned.npz"),
        os.path.join(results_dir, f"heatmap_SV_{dmr_type}_binned.npz"),
        os.path.join(results_dir, f"wholegenome_heatmap_SV_{dmr_type}_binned.npz"),
        os.path.join(results_dir, f"wholegenome_heatmap_{dmr_type}_binned.npz"),
    ]

    for path in candidates:
        if os.path.exists(path):
            return path

    raise ValueError("Heatmap NPZ not found. Checked:\n" + "\n".join(candidates))


# =========================
# Chromosome sizes from SL5 FAI
# =========================
chr_sizes = load_chr_sizes_from_fai(SL5_FAI)
chrom_order = [f"{i:02d}" for i in range(1, 13)]

print("Chromosome sizes loaded from SL5.0.fasta.fai:")
for chrom in chrom_order:
    print(f"{chrom}\t{chr_sizes[chrom]}")


# =========================
# Genome coordinates
# =========================
cumulative_lengths = {}
cum_len = 0
for chrom in chrom_order:
    cumulative_lengths[chrom] = cum_len
    cum_len += chr_sizes[chrom]

genome_length = cum_len
bin_size = 1_000_000
num_bins = genome_length // bin_size + 1


# =========================
# Load SV position map
# =========================
sv_pos_map = build_sv_pos_map(SV_COORD_FILE)
print(f"Loaded {len(sv_pos_map)} SV positions from: {SV_COORD_FILE}")


# =========================
# Build Fisher combined stats
# =========================
if (not PLOT_ONLY) or (not os.path.exists(COMBINED_TSV)):
    input_files = collect_input_files(BASE_DIR, DMR_TYPE, chrom_order)

    if not input_files:
        raise ValueError(f"No input files found for {DMR_TYPE}")

    print(f"Found {len(input_files)} input files")

    keys, key_index_map = build_key_order(input_files[0])
    n_keys = len(keys)

    print(f"Total SV IDs from first file: {n_keys}")
    print("First 10 keys:", keys[:10])

    mapped = sum(1 for k in keys if k in sv_pos_map)
    print(f"SV IDs with coordinate mapping: {mapped}/{n_keys}")

    sum_log_p = np.zeros(n_keys, dtype=np.float64)
    file_count = np.zeros(n_keys, dtype=np.int32)

    for file_idx, file_path in enumerate(input_files, start=1):
        if file_idx % 100 == 0 or file_idx == 1 or file_idx == len(input_files):
            print(f"Processing file {file_idx}/{len(input_files)}: {os.path.basename(file_path)}")

        with gzip.open(file_path, "rt") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 4:
                    continue

                key = parts[0]
                idx = key_index_map.get(key)
                if idx is None:
                    continue

                try:
                    raw_p = float(parts[3]) if parts[3] != "nan" else np.nan
                except ValueError:
                    continue

                p = safe_pvalue(raw_p)
                if p is None:
                    continue

                sum_log_p[idx] += np.log(p)
                file_count[idx] += 1

    combined_chisq = np.zeros(n_keys, dtype=np.float64)
    combined_pvalue = np.ones(n_keys, dtype=np.float64)

    valid = file_count > 0
    combined_chisq[valid] = -2.0 * sum_log_p[valid]
    combined_pvalue[valid] = chi2.sf(combined_chisq[valid], 2 * file_count[valid])

    save_combined_results(
        COMBINED_TSV,
        keys,
        combined_chisq,
        combined_pvalue,
        file_count,
        sv_pos_map,
        cumulative_lengths,
        bin_size
    )

    print(f"Saved combined results: {COMBINED_TSV}")

else:
    print(f"Using existing combined file: {COMBINED_TSV}")


# =========================
# Load combined stats
# =========================
chroms, positions, x_bins, combined_chisq, combined_pvalue, file_count = load_combined_results(COMBINED_TSV)

if len(x_bins) == 0:
    raise ValueError("No combined results loaded. Check SV coordinates and TSV generation.")

manhattan_score = combined_chisq / 1000.0

valid_score = np.isfinite(manhattan_score)
if valid_score.sum() == 0:
    raise ValueError("No valid Manhattan score values")

top1_threshold = np.nanpercentile(manhattan_score[valid_score], 99.0)
top_hits = manhattan_score >= top1_threshold


# =========================
# Load corrected whole-genome heatmap NPZ
# =========================
HEATMAP_NPZ = find_heatmap_npz(RESULTS_DIR, DMR_TYPE)
print(f"Using heatmap NPZ: {HEATMAP_NPZ}")

heatmap_data_binned = load_npz(HEATMAP_NPZ)
print(f"Heatmap NPZ shape: {heatmap_data_binned.shape}, nnz: {heatmap_data_binned.nnz}")

mat = heatmap_data_binned.toarray().astype(np.float32)

plot_mat = np.full(mat.shape, np.nan, dtype=np.float32)
nz = mat > 0
plot_mat[nz] = -np.log10(mat[nz])

if np.all(np.isnan(plot_mat)):
    raise ValueError("Heatmap matrix is empty after filtering")


# =========================
# Axis helpers
# =========================
boundaries = []
mid_bins = []
mid_labels = []

for chrom in chrom_order:
    start_bp = cumulative_lengths[chrom]
    end_bp = start_bp + chr_sizes[chrom]

    start_bin = start_bp // bin_size
    end_bin = end_bp // bin_size
    mid_bin = (start_bin + end_bin) / 2.0

    boundaries.append(start_bin)
    mid_bins.append(mid_bin)
    mid_labels.append(chrom)


# =========================
# Manhattan colors
# =========================
chrom_colors = {
    "01": "#f08070", "02": "#808080",
    "03": "#f08070", "04": "#808080",
    "05": "#f08070", "06": "#808080",
    "07": "#f08070", "08": "#808080",
    "09": "#f08070", "10": "#808080",
    "11": "#f08070", "12": "#808080",
}
point_colors = np.array([chrom_colors.get(c, "#808080") for c in chroms])


# =========================
# Plot final combined figure
# =========================
fig = plt.figure(figsize=(16, 18))
gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1.2, 4.5], hspace=0.05)

ax_top = fig.add_subplot(gs[0, 0])
ax_heat = fig.add_subplot(gs[1, 0], sharex=ax_top)

# ---- top Manhattan-like panel
ax_top.scatter(
    x_bins[~top_hits],
    manhattan_score[~top_hits],
    c=point_colors[~top_hits],
    s=18,
    alpha=0.95,
    linewidths=0
)

ax_top.scatter(
    x_bins[top_hits],
    manhattan_score[top_hits],
    c="red",
    s=28,
    alpha=1.0,
    edgecolors="black",
    linewidths=0.8
)

for b in boundaries:
    ax_top.axvline(
        x=b - 0.5,
        color="lightgrey",
        linestyle="--",
        linewidth=0.6,
        alpha=0.8
    )

ax_top.axhline(
    y=top1_threshold,
    color="blue",
    linestyle="--",
    linewidth=1.2,
    label="Top 1% threshold"
)

ax_top.set_ylabel("Combined chisq (/1000)", fontsize=12)
ax_top.set_title(f"{DMR_TYPE}: combined Fisher SV track + genome heatmap", fontsize=16)
ax_top.legend(loc="upper right", frameon=True)
ax_top.tick_params(axis="x", labelbottom=False)

# ---- bottom heatmap panel
cmap = plt.cm.get_cmap("magma_r").copy()
cmap.set_bad("white")

im = ax_heat.imshow(
    plot_mat,
    origin="lower",
    interpolation="nearest",
    aspect="equal",
    cmap=cmap
)

for b in boundaries:
    ax_heat.axvline(
        x=b - 0.5,
        color="grey",
        linestyle="--",
        linewidth=0.6,
        alpha=0.7
    )
    ax_heat.axhline(
        y=b - 0.5,
        color="grey",
        linestyle="--",
        linewidth=0.6,
        alpha=0.7
    )

ax_heat.set_xticks(mid_bins)
ax_heat.set_xticklabels(mid_labels)
ax_heat.set_yticks(mid_bins)
ax_heat.set_yticklabels(mid_labels)

ax_heat.set_xlabel("SV genome bin (all chromosomes, 1 Mb, SL5)", fontsize=12)
ax_heat.set_ylabel("DMR genome bin (all chromosomes, 1 Mb, SL5)", fontsize=12)

cbar = fig.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.03)
cbar.set_label("-log10(Value)", fontsize=11)

out_png = os.path.join(PLOTS_DIR, f"{DMR_TYPE}_combined_SV_track_plus_heatmap.png")
out_pdf = os.path.join(PLOTS_DIR, f"{DMR_TYPE}_combined_SV_track_plus_heatmap.pdf")

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
