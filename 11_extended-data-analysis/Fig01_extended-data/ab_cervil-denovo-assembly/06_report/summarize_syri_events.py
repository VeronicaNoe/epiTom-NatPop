#!/usr/bin/env python3
from pathlib import Path
import csv
from collections import defaultdict

PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")
INDIR = PROJECT / "06_report/syri_plotsr_chr12"
SYRI = INDIR / "cervil_vs_SL5.syri.out"
FAI = INDIR / "SL5.chr12.renamed.fasta.fai"
OUT = PROJECT / "06_report/syri_event_summary.tsv"

MIN_LEN = 50000  # change to 100000 for stricter filtering

event_keywords = {
    "SYN", "INV", "TRANS", "DUP", "INVDP", "INVTR", "NOTAL", "SNP", "INS", "DEL", "HDR"
}

ref_sizes = {}
with open(FAI) as f:
    for line in f:
        p = line.rstrip().split("\t")
        ref_sizes[p[0]] = int(p[1])

total_ref = sum(ref_sizes.values())

events = []

with open(SYRI) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip().split("\t")

        if len(p) < 4:
            continue

        ref_chr = p[0]
        try:
            ref_start = int(p[1])
            ref_end = int(p[2])
        except ValueError:
            continue

        length = abs(ref_end - ref_start) + 1

        # SyRI event type can vary by output version; scan fields for known labels
        evtype = None
        for field in reversed(p):
            if field in event_keywords:
                evtype = field
                break

        if evtype is None:
            continue

        if length < MIN_LEN:
            continue

        events.append((ref_chr, min(ref_start, ref_end), max(ref_start, ref_end), evtype, length))

# merge intervals per event type
by_type = defaultdict(list)
for chrom, start, end, evtype, length in events:
    by_type[evtype].append((chrom, start, end))

rows = []
for evtype, intervals in by_type.items():
    intervals.sort()
    merged = []

    for chrom, start, end in intervals:
        if not merged or merged[-1][0] != chrom or start > merged[-1][2] + 1:
            merged.append([chrom, start, end])
        else:
            merged[-1][2] = max(merged[-1][2], end)

    bp = sum(end - start + 1 for chrom, start, end in merged)
    rows.append({
        "event_type": evtype,
        "n_raw_events_ge_50kb": len(intervals),
        "n_merged_intervals": len(merged),
        "ref_bp_covered": bp,
        "ref_pct_covered": 100 * bp / total_ref,
    })

rows.sort(key=lambda x: x["ref_bp_covered"], reverse=True)

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=["event_type", "n_raw_events_ge_50kb", "n_merged_intervals", "ref_bp_covered", "ref_pct_covered"],
        delimiter="\t"
    )
    w.writeheader()
    for r in rows:
        w.writerow(r)

print(f"Wrote: {OUT}")
