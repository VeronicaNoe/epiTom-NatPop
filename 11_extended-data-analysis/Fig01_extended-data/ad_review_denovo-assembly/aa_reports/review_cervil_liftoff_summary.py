#!/usr/bin/env python3
from pathlib import Path
import csv

PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")
GFF = PROJECT / "05_annotation/liftoff_final/cervil.final.liftoff.gff3"
UNMAPPED = PROJECT / "05_annotation/liftoff_final/cervil.final.liftoff.unmapped.txt"
OUT = PROJECT / "06_report/liftoff_summary.tsv"

total = 0
genes = 0
mrna = 0
cds = 0

with open(GFF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip().split("\t")
        if len(parts) < 9:
            continue
        total += 1
        feat = parts[2]
        if feat == "gene":
            genes += 1
        elif feat in ("mRNA", "transcript"):
            mrna += 1
        elif feat == "CDS":
            cds += 1

unmapped = 0
if UNMAPPED.exists():
    with open(UNMAPPED) as f:
        for line in f:
            if line.strip():
                unmapped += 1

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerow(["gff_features_total", total])
    w.writerow(["genes_mapped", genes])
    w.writerow(["transcripts_mapped", mrna])
    w.writerow(["cds_features", cds])
    w.writerow(["genes_unmapped", unmapped])

print(f"Wrote: {OUT}")
