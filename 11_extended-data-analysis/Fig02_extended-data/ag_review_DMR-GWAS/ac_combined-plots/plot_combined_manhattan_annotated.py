#!/users/bioinfo/vibanez/anaconda3/bin/python

import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
from scipy.sparse import load_npz, save_npz


# =========================
# Arguments
# =========================
parser = argparse.ArgumentParser()
parser.add_argument("dmr_type", choices=["C-DMR", "CG-DMR"])
parser.add_argument("--bin-size", type=int, default=1_000_000)
parser.add_argument("--top-frac", type=float, default=0.01)
parser.add_argument("--max-labels", type=int, default=40)
args = parser.parse_args()

DMR_TYPE = args.dmr_type
BIN_SIZE = args.bin_size
TOP_FRAC = args.top_frac
MAX_LABELS = args.max_labels


# =========================
# Embedded paths
# =========================
SL25_FAI = "/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/ITAG2.4_genomic.fasta.fai"

# SNP
SNP_TABLE = f"/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/ab_DMR-GWAS-SNPs/combined_fisher_SNP_{DMR_TYPE}.tsv"
SNP_NPZ   = f"/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/ab_DMR-GWAS-SNPs/heatmap_{DMR_TYPE}_SNP_SL25_binned.npz"

# TIP
TIP_TABLE = f"/mnt/disk2/vibanez/otherAnalysis/11_DMR-GWAS-TIPs/ae_data-analysis/results/combined_fisher_{DMR_TYPE}.tsv"
TIP_NPZ   = f"/mnt/disk2/vibanez/otherAnalysis/11_DMR-GWAS-TIPs/ae_data-analysis/results/heatmap_{DMR_TYPE}_binned.npz"

# SV
SV_TABLE = f"/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/aa_DMR-GWAS-SV_SL25/combined_fisher_SV_{DMR_TYPE}.tsv"
SV_NPZ   = f"/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/aa_DMR-GWAS-SV_SL25/heatmap_{DMR_TYPE}_SV_SL25_binned.npz"

# SNP annotation table
SNP_ANNOT = "/mnt/disk2/vibanez/10_data-analysis/Fig3/ab_data-analysis/results/03.9_SNPs_over_epiGenes_annotated.tsv"

OUTDIR = "/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/ac_combined-plots"
os.makedirs(OUTDIR, exist_ok=True)


# =========================
# Scale raw chisq for display / diagnostics
# =========================
SCORE_DIVISOR = {
    "SNP": 1000.0,
    "TIP": 1000.0,
    "SV": 1000.0,
}


# =========================
# Helpers
# =========================
PRIMARY_CHROMS = [f"{i:02d}" for i in range(1, 13)]

def normalize_chr(chrn):
    s = str(chrn).strip()
    s = re.sub(r"^SL2\.50ch", "", s, flags=re.IGNORECASE)
    s = re.sub(r"^(chr|ch)", "", s, flags=re.IGNORECASE)
    m = re.search(r"(\d+)$", s)
    if m is None:
        return None
    return f"{int(m.group(1)):02d}"

def chrom_alt_color(chrom):
    return "#808080" if (int(chrom) % 2 == 0) else "#f08070"

def load_chr_sizes_from_fai(fai_file):
    chr_sizes = {}
    with open(fai_file, "rt") as f:
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
        raise ValueError(f"Missing chr sizes in FAI: {missing}")

    return chr_sizes

def parse_key_to_chr_pos(key):
    key = str(key)
    if ":" in key:
        parts = key.split(":")
        if len(parts) >= 2:
            chrom = normalize_chr(parts[0])
            try:
                pos = int(parts[1])
            except ValueError:
                pos = None
            return chrom, pos
    return None, None

def natural_key(s):
    return [
        int(x) if x.isdigit() else x
        for x in re.split(r"(\d+)", str(s))
    ]

def clean_text(x):
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if s.lower() in {"", "nan", "none", "na"}:
        return ""
    return s

def clean_gene_label(x):
    """
    Clean the gene symbol column.

    Examples:
      DRM2,DRM1       -> DRM1/2
      AGO2,AGO1       -> AGO1/2
      NRPE9A,NRPE9B   -> NRPE9A/B
      DMS11/MORC6     -> DMS11/MORC6
    """
    s = clean_text(x)
    if s == "":
        return np.nan

    # Only split comma-separated alternative names.
    # Keep existing slash labels such as DMS11/MORC6.
    tokens = [t.strip() for t in re.split(r"\s*,\s*", s) if t.strip()]
    if len(tokens) == 0:
        return np.nan
    if len(tokens) == 1:
        return tokens[0]

    tokens = sorted(set(tokens), key=natural_key)

    common = os.path.commonprefix(tokens)

    # Compress only if the shared prefix is meaningful.
    if len(common) >= 2 and all(len(t) > len(common) for t in tokens):
        suffixes = [t[len(common):] for t in tokens]
        return common + "/".join(suffixes)

    return "/".join(tokens)

def collapse_unique_labels(vals):
    labels = []
    seen = set()
    for v in vals:
        if pd.isna(v):
            continue
        s = str(v).strip()
        if s == "":
            continue
        if s not in seen:
            labels.append(s)
            seen.add(s)

    if len(labels) == 0:
        return np.nan
    return ";".join(labels)

def read_snp_annotation_table(path):
    annot = pd.read_csv(path, sep="\t")

    required = {"SNPs", "Name"}
    missing = required - set(annot.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")

    chr_pos = annot["SNPs"].apply(parse_key_to_chr_pos)
    annot["chrom"] = [x[0] for x in chr_pos]
    annot["pos"]   = [x[1] for x in chr_pos]

    annot = annot.dropna(subset=["chrom", "pos"]).copy()
    annot["chrom"] = annot["chrom"].apply(normalize_chr)
    annot["pos"] = annot["pos"].astype(int)

    annot["gene_label"] = annot["Name"].apply(clean_gene_label)

    annot = annot.dropna(subset=["gene_label"]).copy()

    annot = (
        annot
        .groupby(["chrom", "pos"], as_index=False)
        .agg(gene_label=("gene_label", collapse_unique_labels))
    )

    return annot

def read_combined_table(path, marker_type, cumulative_lengths, bin_size):
    dt = pd.read_csv(path, sep="\t")

    if "chrom" not in dt.columns or "pos" not in dt.columns:
        if "key" not in dt.columns:
            raise ValueError(f"{path} has neither chrom/pos nor key columns")
        chr_pos = dt["key"].apply(parse_key_to_chr_pos)
        dt["chrom"] = [x[0] for x in chr_pos]
        dt["pos"]   = [x[1] for x in chr_pos]

    needed = {"key", "combined_chisq", "combined_pvalue", "file_count", "chrom", "pos"}
    missing = needed - set(dt.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")

    dt = dt.copy()
    dt["chrom"] = dt["chrom"].apply(normalize_chr)
    dt = dt[dt["chrom"].isin(PRIMARY_CHROMS)].copy()

    dt["pos"] = pd.to_numeric(dt["pos"], errors="coerce")
    dt = dt.dropna(subset=["pos"])
    dt["pos"] = dt["pos"].astype(int)

    dt["combined_chisq"] = pd.to_numeric(dt["combined_chisq"], errors="coerce")
    dt["combined_pvalue"] = pd.to_numeric(dt["combined_pvalue"], errors="coerce")
    dt["file_count"] = pd.to_numeric(dt["file_count"], errors="coerce")

    dt = dt.replace([np.inf, -np.inf], np.nan)
    dt = dt.dropna(subset=["combined_chisq", "combined_pvalue", "file_count"])

    dt["marker_type"] = marker_type
    dt["linear_pos"] = [cumulative_lengths[c] + p for c, p in zip(dt["chrom"], dt["pos"])]
    dt["x_bin"] = dt["linear_pos"] / float(bin_size)

    divisor = SCORE_DIVISOR[marker_type]
    dt["raw_score"] = dt["combined_chisq"] / divisor
    dt["chr_color"] = dt["chrom"].apply(chrom_alt_color)

    return dt

def add_within_class_tail_score(dt):
    dt = dt.copy()

    dt["rank_within_type"] = (
        dt.groupby("marker_type")["raw_score"]
        .rank(method="average", ascending=False)
    )

    dt["n_within_type"] = (
        dt.groupby("marker_type")["raw_score"]
        .transform("size")
    )

    dt["tail_fraction_within_type"] = (
        dt["rank_within_type"] / (dt["n_within_type"] + 1.0)
    )

    dt["plot_score"] = -np.log10(dt["tail_fraction_within_type"])

    return dt

def add_snp_annotations(merged_tab, annot_path):
    annot = read_snp_annotation_table(annot_path)

    out = merged_tab.merge(
        annot,
        on=["chrom", "pos"],
        how="left"
    )

    # Only SNPs should keep SNP gene annotations.
    out.loc[out["marker_type"] != "SNP", "gene_label"] = np.nan

    n_annot_total = len(annot)
    n_annot_matched = out.loc[
        (out["marker_type"] == "SNP") & out["gene_label"].notna()
    ].shape[0]

    print(f"SNP annotation rows after cleaning: {n_annot_total}")
    print(f"Annotated SNPs matched in Manhattan table: {n_annot_matched}")

    if n_annot_matched == 0:
        print("WARNING: no annotated SNPs matched. Check coordinate system or SNP key format.")

    return out

def merge_npzs(npz_paths):
    mats = [load_npz(p).tocsr() for p in npz_paths]
    shape0 = mats[0].shape

    for i, m in enumerate(mats[1:], start=2):
        if m.shape != shape0:
            raise ValueError(f"NPZ shape mismatch: matrix1={shape0}, matrix{i}={m.shape}")

    merged = mats[0].tolil()

    for m in mats[1:]:
        coo = m.tocoo()
        for r, c, v in zip(coo.row, coo.col, coo.data):
            cur = merged[r, c]
            if cur == 0 or v < cur:
                merged[r, c] = v

    return merged.tocsr()

def annotate_snp_labels(ax, annot_dt, y_upper):
    """
    Add staggered SNP labels with arrows.
    """
    if len(annot_dt) == 0:
        return

    annot_dt = annot_dt.sort_values("x_bin").reset_index(drop=True)

    dx_cycle = [-20, -10, 0, 10, 20]
    dy_cycle = [18, 28, 38, 48, 58]

    for i, row in annot_dt.iterrows():
        dx = dx_cycle[i % len(dx_cycle)]
        dy = dy_cycle[i % len(dy_cycle)]

        ax.annotate(
            row["gene_label"],
            xy=(row["x_bin"], row["plot_score"]),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            fontstyle="italic",
            arrowprops=dict(
                arrowstyle="-",
                color="black",
                linewidth=0.5,
                shrinkA=0,
                shrinkB=0
            ),
            annotation_clip=False,
            clip_on=False
        )


# =========================
# Check files
# =========================
for p in [SL25_FAI, SNP_TABLE, SNP_NPZ, TIP_TABLE, TIP_NPZ, SV_TABLE, SV_NPZ, SNP_ANNOT]:
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing file: {p}")


# =========================
# Genome layout
# =========================
chr_sizes = load_chr_sizes_from_fai(SL25_FAI)
chrom_order = [f"{i:02d}" for i in range(1, 13)]

cumulative_lengths = {}
cum = 0
for chrom in chrom_order:
    cumulative_lengths[chrom] = cum
    cum += chr_sizes[chrom]

boundaries = []
mid_bins = []
mid_labels = []

for chrom in chrom_order:
    start_bp = cumulative_lengths[chrom]
    end_bp = start_bp + chr_sizes[chrom]

    start_bin = start_bp // BIN_SIZE
    end_bin = end_bp // BIN_SIZE
    mid_bin = (start_bin + end_bin) / 2.0

    boundaries.append(start_bin)
    mid_bins.append(mid_bin)
    mid_labels.append(str(int(chrom)))


# =========================
# Merge combined tables
# =========================
snp = read_combined_table(SNP_TABLE, "SNP", cumulative_lengths, BIN_SIZE)
tip = read_combined_table(TIP_TABLE, "TIP", cumulative_lengths, BIN_SIZE)
sv  = read_combined_table(SV_TABLE,  "SV",  cumulative_lengths, BIN_SIZE)

merged_tab = pd.concat([snp, tip, sv], ignore_index=True)
merged_tab = add_within_class_tail_score(merged_tab)
merged_tab = add_snp_annotations(merged_tab, SNP_ANNOT)

merged_tab.to_csv(
    os.path.join(OUTDIR, f"merged_marker_tables_{DMR_TYPE}_SL25.with_SNP_annotations.tsv"),
    sep="\t",
    index=False
)

print("Loaded rows:")
print(f"  SNP: {len(snp)}")
print(f"  TIP: {len(tip)}")
print(f"  SV : {len(sv)}")

for marker_type in ["SNP", "TIP", "SV"]:
    sub = merged_tab[merged_tab["marker_type"] == marker_type]
    if len(sub) > 0:
        print(
            f"{marker_type} raw_score range: "
            f"min={sub['raw_score'].min():.4f}, "
            f"median={sub['raw_score'].median():.4f}, "
            f"max={sub['raw_score'].max():.4f}"
        )
        print(
            f"{marker_type} plot_score range: "
            f"min={sub['plot_score'].min():.4f}, "
            f"median={sub['plot_score'].median():.4f}, "
            f"max={sub['plot_score'].max():.4f}"
        )


# =========================
# Manhattan
# =========================
marker_shapes = {
    "SNP": "o",
    "TIP": "D",
    "SV":  "s"
}

marker_sizes = {
    "SNP": 10,
    "TIP": 16,
    "SV":  16
}

fill_alpha = {
    "SNP": 0.35,
    "TIP": 0.55,
    "SV":  0.55
}

top_line_y = -np.log10(TOP_FRAC)

fig, ax = plt.subplots(figsize=(16, 5))

for marker_type in ["SNP", "TIP", "SV"]:
    sub = merged_tab[merged_tab["marker_type"] == marker_type].copy()
    facecols = [to_rgba(c, alpha=fill_alpha[marker_type]) for c in sub["chr_color"]]

    ax.scatter(
        sub["x_bin"],
        sub["plot_score"],
        facecolors=facecols,
        edgecolors="black",
        linewidths=0.15,
        marker=marker_shapes[marker_type],
        s=marker_sizes[marker_type],
        rasterized=True
    )

# Top 1% horizontal line
ax.axhline(
    y=top_line_y,
    color="blue",
    linestyle=":",
    linewidth=1.0,
    alpha=0.9
)

xmax = merged_tab["x_bin"].max()
ax.text(
    xmax,
    top_line_y + 0.03,
    f"top {int(TOP_FRAC * 100)}%",
    ha="right",
    va="bottom",
    fontsize=9,
    color="black"
)

# Chromosome boundaries
for b in boundaries:
    ax.axvline(
        x=b - 0.5,
        color="lightgrey",
        linestyle="--",
        linewidth=0.6,
        alpha=0.8
    )

# Select all annotated SNPs
annot_dt = merged_tab[
    (merged_tab["marker_type"] == "SNP") &
    (merged_tab["gene_label"].notna())
].copy()

# Avoid duplicate labels if the same SNP appears more than once
annot_dt = annot_dt.drop_duplicates(subset=["chrom", "pos", "gene_label"])

# Plot labels in genomic order
annot_dt = annot_dt.sort_values(["chrom", "pos"]).reset_index(drop=True)

print(f"Annotated SNPs plotted/highlighted: {len(annot_dt)}")

# Highlight annotated SNPs
if len(annot_dt) > 0:
    ax.scatter(
        annot_dt["x_bin"],
        annot_dt["plot_score"],
        facecolors="none",
        edgecolors="black",
        linewidths=0.8,
        marker="o",
        s=45,
        zorder=5
    )

ymax = merged_tab["plot_score"].max()
y_upper = max(ymax * 1.25, top_line_y * 1.25)
ax.set_ylim(0, y_upper)

annotate_snp_labels(ax, annot_dt, y_upper)

ax.set_xticks(mid_bins)
ax.set_xticklabels(mid_labels, rotation=0)
ax.set_xlabel("Chromosome")
ax.set_ylabel("-log10(rank / (N+1)) within marker type")
ax.set_title(f"{DMR_TYPE}: merged SNP + TIP + SV Manhattan with SNP gene annotations")

legend_handles = [
    Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor="black",
           linestyle="None", markersize=7, label="SNP"),
    Line2D([0], [0], marker="D", markerfacecolor="white", markeredgecolor="black",
           linestyle="None", markersize=7, label="TIP"),
    Line2D([0], [0], marker="s", markerfacecolor="white", markeredgecolor="black",
           linestyle="None", markersize=7, label="SV"),
]
ax.legend(handles=legend_handles, frameon=True, loc="upper right")

plt.tight_layout()
plt.savefig(
    os.path.join(OUTDIR, f"{DMR_TYPE}_merged_SNP_TIP_SV_manhattan_SL25.annotated.png"),
    dpi=300
)
plt.savefig(
    os.path.join(OUTDIR, f"{DMR_TYPE}_merged_SNP_TIP_SV_manhattan_SL25.annotated.pdf"),
    dpi=300
)
plt.close()


# =========================
# Merge NPZ heatmaps
# =========================
merged_npz = merge_npzs([SNP_NPZ, TIP_NPZ, SV_NPZ])

merged_npz_path = os.path.join(
    OUTDIR,
    f"heatmap_{DMR_TYPE}_SNP_TIP_SV_SL25_binned.npz"
)
save_npz(merged_npz_path, merged_npz)

print(f"Saved merged NPZ: {merged_npz_path}")


# =========================
# Heatmap
# =========================
mat = merged_npz.toarray().astype(np.float32)
plot_mat = np.full(mat.shape, np.nan, dtype=np.float32)
nz = mat > 0
plot_mat[nz] = -np.log10(mat[nz])

cmap = plt.cm.get_cmap("magma_r").copy()
cmap.set_bad("white")

fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(
    plot_mat,
    origin="lower",
    interpolation="nearest",
    aspect="equal",
    cmap=cmap
)

for b in boundaries:
    ax.axvline(x=b - 0.5, color="grey", linestyle="--", linewidth=0.6, alpha=0.7)
    ax.axhline(y=b - 0.5, color="grey", linestyle="--", linewidth=0.6, alpha=0.7)

ax.set_xticks(mid_bins)
ax.set_xticklabels(mid_labels)
ax.set_yticks(mid_bins)
ax.set_yticklabels(mid_labels)

ax.set_xlabel("Marker genome bin (SL2.5, 1 Mb)")
ax.set_ylabel("DMR genome bin (SL2.5, 1 Mb)")
ax.set_title(f"{DMR_TYPE}: merged SNP + TIP + SV heatmap (SL2.5)")

cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.03)
cbar.set_label("-log10(min p)")

plt.tight_layout()
plt.savefig(
    os.path.join(OUTDIR, f"{DMR_TYPE}_merged_SNP_TIP_SV_heatmap_SL25.png"),
    dpi=300
)
plt.savefig(
    os.path.join(OUTDIR, f"{DMR_TYPE}_merged_SNP_TIP_SV_heatmap_SL25.pdf"),
    dpi=300
)
plt.close()

