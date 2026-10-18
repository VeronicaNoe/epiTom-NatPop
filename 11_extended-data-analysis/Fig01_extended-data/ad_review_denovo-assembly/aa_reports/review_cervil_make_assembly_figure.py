from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")
REPORT = PROJECT / "06_report"

SUMMARY = REPORT / "assembly_summary.tsv"
PAF = REPORT / "ragtag.vs_SL5.paf"
RAGTAG_STATS = PROJECT / "03_ragtag/ragtag.scaffold.stats"

OUT = REPORT / "cervil_assembly_overview.pdf"

df = pd.read_csv(SUMMARY, sep="\t")

# -------------------------
# Parse RagTag stats
# -------------------------
rag = {}
with open(RAGTAG_STATS) as f:
    header = f.readline().strip().split("\t")
    values = f.readline().strip().split("\t")
    rag = dict(zip(header, map(float, values)))

placed_bp = rag["placed_bp"]
unplaced_bp = rag["unplaced_bp"]
gap_bp = rag["gap_bp"]

# -------------------------
# Parse PAF for dotplot
# -------------------------
x, y = [], []
with open(PAF) as f:
    for line in f:
        p = line.rstrip().split("\t")
        if len(p) < 12:
            continue

        qstart, qend = int(p[2]), int(p[3])
        tstart, tend = int(p[7]), int(p[8])

        x.extend([tstart, tend, None])
        y.extend([qstart, qend, None])

# -------------------------
# Figure
# -------------------------
fig = plt.figure(figsize=(14, 10))

ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

# -------------------------
# A) Assembly metrics
# -------------------------
plot_df = df.copy()
plot_df["assembly_label"] = plot_df["assembly"].replace({
    "hifiasm_raw": "hifiasm",
    "dorado_polish": "polished",
    "ragtag": "RagTag"
})

ax1.bar(plot_df["assembly_label"], plot_df["n50_bp"] / 1e6)
ax1.set_ylabel("N50 (Mb)")
ax1.set_title("A) Assembly contiguity")

for i, row in plot_df.iterrows():
    ax1.text(
        i,
        row["n50_bp"] / 1e6,
        f"L50={int(row['l50'])}\n{int(row['n_seqs'])} seqs",
        ha="center",
        va="bottom",
        fontsize=9,
    )

# -------------------------
# B) BUSCO stacked bars
# -------------------------
busco_cols = [
    "busco_single_pct",
    "busco_duplicated_pct",
    "busco_fragmented_pct",
    "busco_missing_pct",
]

bottom = [0] * len(plot_df)
labels = ["Single-copy", "Duplicated", "Fragmented", "Missing"]

for col, lab in zip(busco_cols, labels):
    vals = plot_df[col].astype(float).values
    ax2.barh(plot_df["assembly_label"], vals, left=bottom, label=lab)
    bottom = [b + v for b, v in zip(bottom, vals)]

ax2.set_xlim(0, 100)
ax2.set_xlabel("BUSCO genes (%)")
ax2.set_title("B) BUSCO completeness")
ax2.legend(frameon=False, fontsize=9, loc="lower right")

# -------------------------
# C) Synteny dotplot
# -------------------------
ax3.plot(x, y, linewidth=0.2)
ax3.set_xlabel("SL5 coordinate")
ax3.set_ylabel("Cervil RagTag coordinate")
ax3.set_title("C) Cervil–SL5 whole-genome alignment")

# -------------------------
# D) RagTag anchoring
# -------------------------
vals = [placed_bp / 1e6, unplaced_bp / 1e6, gap_bp / 1e6]
labs = ["Anchored", "Unplaced", "Gaps"]

ax4.bar(labs, vals)
ax4.set_ylabel("Mb")
ax4.set_title("D) RagTag anchoring summary")

total = placed_bp + unplaced_bp + gap_bp
for i, v in enumerate(vals):
    pct = 100 * (v * 1e6) / total
    ax4.text(i, v, f"{pct:.1f}%", ha="center", va="bottom", fontsize=10)

plt.tight_layout()
plt.savefig(OUT)
plt.close()

print(f"Wrote: {OUT}")
