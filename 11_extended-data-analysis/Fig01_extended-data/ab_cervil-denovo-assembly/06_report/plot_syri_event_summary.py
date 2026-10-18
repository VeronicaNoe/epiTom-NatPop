#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")
REPORT = PROJECT / "06_report"

IN = REPORT / "syri_event_summary.filtered_100kb.tsv"
OUT = REPORT / "syri_event_summary.filtered_100kb.pdf"

df = pd.read_csv(IN, sep="\t")

# Convert bp to Mb
df["merged_ref_Mb"] = df["merged_ref_bp"] / 1e6

# Order by bp covered
df = df.sort_values("merged_ref_bp", ascending=True)

labels = [
    f"{row.event_type}\n{row.confidence_class}"
    for row in df.itertuples()
]

plt.figure(figsize=(8, 5))
plt.barh(labels, df["merged_ref_Mb"])
plt.xlabel("SL5 reference sequence covered (Mb)")
plt.ylabel("SyRI class")
plt.title("Large SyRI events (≥100 kb)")
plt.tight_layout()
plt.savefig(OUT)
plt.close()

print(f"Wrote: {OUT}")
