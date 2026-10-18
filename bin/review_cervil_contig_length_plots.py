#!/usr/bin/env python3

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")
REPORT = PROJECT / "06_report"

ASSEMBLIES = {
    "hifiasm": PROJECT / "01_hifiasm/cervil.primary.fasta",
    "polished": PROJECT / "02_polish/cervil.polished.fasta",
    "RagTag": PROJECT / "03_ragtag/ragtag.scaffold.fasta",
}

OUT = REPORT / "cervil_contig_length_distribution.pdf"

def fasta_lengths(path):
    lengths = []
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if n > 0:
                    lengths.append(n)
                n = 0
            else:
                n += len(line.strip())
        if n > 0:
            lengths.append(n)
    return np.array(lengths)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# A) Histogram
for label, fasta in ASSEMBLIES.items():
    lengths_bp = fasta_lengths(fasta)
    lengths_mb = lengths_bp / 1e6

    axes[0].hist(
        lengths_mb,
        bins=np.logspace(np.log10(0.001), np.log10(100), 70),
        histtype="step",
        linewidth=1.5,
        label=label,
    )

axes[0].axvline(0.001, linestyle="--", linewidth=0.8)
axes[0].axvline(1, linestyle="--", linewidth=0.8)
axes[0].set_xscale("log")
axes[0].set_xlim(0.001, 100)
axes[0].set_xlabel("Sequence length (Mb; 0.001 Mb = 1 kb, 1 Mb marked)")
axes[0].set_ylabel("Number of sequences")
axes[0].set_title("A) Sequence length distribution")
axes[0].legend(frameon=False)

# B) Cumulative curve
for label, fasta in ASSEMBLIES.items():
    lengths = np.sort(fasta_lengths(fasta))[::-1]
    cum = np.cumsum(lengths) / 1e6
    x = np.arange(1, len(lengths) + 1)

    total = cum[-1]
    half = total * 0.5
    n50_idx = np.searchsorted(cum, half)
    n50_rank = n50_idx + 1

    axes[1].plot(x, cum, linewidth=1.5, label=label)
    axes[1].axvline(n50_rank, linestyle="--", linewidth=0.8)

axes[1].set_xscale("log")
axes[1].set_xlabel("Top N longest sequences, ranked by length")
axes[1].set_ylabel("Cumulative assembly size (Mb)")
axes[1].set_title("B) Cumulative assembly curve")
axes[1].legend(frameon=False)

plt.tight_layout()
plt.savefig(OUT)
plt.close()

print(f"Wrote: {OUT}")
