#!/usr/bin/env python3

from pathlib import Path
import csv
from collections import defaultdict

PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")

SYRI = PROJECT / "06_report/syri_plotsr_chr12/cervil_vs_SL5.syri.out"
REF_FAI = PROJECT / "06_report/syri_plotsr_chr12/SL5.chr12.renamed.fasta.fai"
AGP = PROJECT / "03_ragtag/ragtag.scaffold.agp"

# This is the EDTA file you showed.
# Note: this is intact TE annotation, not necessarily all repeats.
EDTA_GFF = PROJECT / "05_annotation/edta_final/ragtag.scaffold.fasta.mod.EDTA.intact.gff3"

OUT_EVENTS = PROJECT / "06_report/syri_events.filtered_100kb.tsv"
OUT_SUMMARY = PROJECT / "06_report/syri_event_summary.filtered_100kb.tsv"

MIN_LEN = 100_000
JOIN_WINDOW = 10_000
TE_FRAC_CUTOFF = 0.50

EVENT_TYPES = {"SYN", "INV", "TRANS", "DUP", "INVDP", "INVTR", "HDR", "NOTAL", "DEL", "INS"}


def chr_to_ragtag(chrom):
    # Chr10 -> 10_RagTag
    if chrom.startswith("Chr"):
        return chrom.replace("Chr", "") + "_RagTag"
    return chrom


def ragtag_to_chr(chrom):
    # 10_RagTag -> Chr10
    if chrom.endswith("_RagTag"):
        return "Chr" + chrom.replace("_RagTag", "")
    return chrom


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for chrom, start, end in intervals[1:]:
        last = merged[-1]
        if chrom == last[0] and start <= last[2] + 1:
            last[2] = max(last[2], end)
        else:
            merged.append([chrom, start, end])
    return merged


def overlap_bp(interval, intervals):
    chrom, start, end = interval
    ov = 0
    for c, s, e in intervals:
        if c != chrom:
            continue
        if e < start:
            continue
        if s > end:
            break
        ov += max(0, min(end, e) - max(start, s) + 1)
    return ov


# Reference genome size
ref_sizes = {}
with open(REF_FAI) as f:
    for line in f:
        p = line.rstrip().split("\t")
        ref_sizes[p[0]] = int(p[1])
total_ref_bp = sum(ref_sizes.values())


# RagTag gap/join positions from AGP
joins = defaultdict(list)

with open(AGP) as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip().split("\t")
        if len(p) < 9:
            continue

        obj = p[0]                 # e.g. 10_RagTag
        obj_chr = ragtag_to_chr(obj)
        start = int(p[1])
        end = int(p[2])
        comp_type = p[4]

        if comp_type == "U":
            joins[obj_chr].append((obj_chr, start, end))

for chrom in joins:
    joins[chrom] = merge_intervals(joins[chrom])


# EDTA repeat intervals, converted to ChrN names
te_intervals = defaultdict(list)

if EDTA_GFF.exists():
    with open(EDTA_GFF) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip().split("\t")
            if len(p) < 5:
                continue
            chrom = ragtag_to_chr(p[0])
            try:
                start = int(p[3])
                end = int(p[4])
            except ValueError:
                continue
            te_intervals[chrom].append((chrom, min(start, end), max(start, end)))

for chrom in te_intervals:
    te_intervals[chrom] = merge_intervals(te_intervals[chrom])


events = []

with open(SYRI) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        p = line.rstrip().split("\t")
        if len(p) < 11:
            continue

        ref_chr = p[0]
        try:
            ref_start = int(p[1])
            ref_end = int(p[2])
        except ValueError:
            continue

        ref_start, ref_end = sorted((ref_start, ref_end))
        ref_len = ref_end - ref_start + 1

        evtype = p[10]
        if evtype not in EVENT_TYPES:
            continue

        if ref_len < MIN_LEN:
            continue

        qry_chr = p[5]
        q_start = p[6]
        q_end = p[7]

        has_query = qry_chr != "-" and q_start != "-" and q_end != "-"

        te_bp = 0
        te_frac = 0.0
        near_join = False
        qry_len = ""

        if has_query:
            try:
                q_start = int(q_start)
                q_end = int(q_end)
                q_start, q_end = sorted((q_start, q_end))
                qry_len = q_end - q_start + 1

                te_bp = overlap_bp((qry_chr, q_start, q_end), te_intervals.get(qry_chr, []))
                te_frac = te_bp / qry_len if qry_len > 0 else 0.0

                # join proximity: expand event interval by JOIN_WINDOW
                join_ov = overlap_bp(
                    (qry_chr, max(1, q_start - JOIN_WINDOW), q_end + JOIN_WINDOW),
                    joins.get(qry_chr, [])
                )
                near_join = join_ov > 0

            except ValueError:
                has_query = False

        if evtype == "SYN":
            confidence = "syntenic"
        elif not has_query:
            confidence = "no_query_interval"
        elif near_join and te_frac >= TE_FRAC_CUTOFF:
            confidence = "low_confidence_repeat_and_join"
        elif near_join:
            confidence = "near_scaffold_join"
        elif te_frac >= TE_FRAC_CUTOFF:
            confidence = "repeat_associated"
        else:
            confidence = "high_confidence_candidate"

        events.append({
            "ref_chr": ref_chr,
            "ref_start": ref_start,
            "ref_end": ref_end,
            "ref_len_bp": ref_len,
            "qry_chr": qry_chr,
            "qry_start": q_start if has_query else "",
            "qry_end": q_end if has_query else "",
            "qry_len_bp": qry_len,
            "event_type": evtype,
            "te_bp_on_query": te_bp,
            "te_frac_on_query": round(te_frac, 4),
            "near_ragtag_join_10kb": near_join,
            "confidence_class": confidence,
        })


# Write event table
OUT_EVENTS.parent.mkdir(parents=True, exist_ok=True)

event_fields = [
    "ref_chr", "ref_start", "ref_end", "ref_len_bp",
    "qry_chr", "qry_start", "qry_end", "qry_len_bp",
    "event_type",
    "te_bp_on_query", "te_frac_on_query",
    "near_ragtag_join_10kb",
    "confidence_class",
]

with open(OUT_EVENTS, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=event_fields, delimiter="\t")
    w.writeheader()
    for e in events:
        w.writerow(e)


# Summary: merged reference bp by event type and confidence class
summary = defaultdict(list)

for e in events:
    key = (e["event_type"], e["confidence_class"])
    summary[key].append((e["ref_chr"], e["ref_start"], e["ref_end"]))

summary_rows = []

for (evtype, conf), intervals in summary.items():
    merged = merge_intervals(intervals)
    bp = sum(e - s + 1 for c, s, e in merged)
    summary_rows.append({
        "event_type": evtype,
        "confidence_class": conf,
        "n_raw_events_ge_100kb": len(intervals),
        "n_merged_intervals": len(merged),
        "merged_ref_bp": bp,
        "merged_ref_pct": round(100 * bp / total_ref_bp, 4),
    })

summary_rows.sort(key=lambda x: (x["event_type"], x["confidence_class"]))

summary_fields = [
    "event_type", "confidence_class",
    "n_raw_events_ge_100kb",
    "n_merged_intervals",
    "merged_ref_bp",
    "merged_ref_pct",
]

with open(OUT_SUMMARY, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=summary_fields, delimiter="\t")
    w.writeheader()
    for r in summary_rows:
        w.writerow(r)

print(f"Wrote: {OUT_EVENTS}")
print(f"Wrote: {OUT_SUMMARY}")
