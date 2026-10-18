#!/usr/bin/env python3
from pathlib import Path
import csv

PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")
EDTA_DIR = PROJECT / "05_annotation/edta_final"
OUT = PROJECT / "06_report/edta_summary.tsv"

sum_files = list(EDTA_DIR.glob("*.EDTA.TEanno.sum")) + list(EDTA_DIR.glob("*.mod.EDTA.TEanno.sum"))

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="") as fout:
    w = csv.writer(fout, delimiter="\t")
    w.writerow(["source_file", "line"])
    for sf in sum_files:
        with open(sf) as f:
            for line in f:
                line = line.rstrip()
                if line:
                    w.writerow([sf.name, line])

print(f"Wrote: {OUT}")
