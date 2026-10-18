#!/usr/bin/env python3
from pathlib import Path
import csv

PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")
REPORT = PROJECT / "06_report"

EVENTS = REPORT / "syri_events.filtered_100kb.tsv"
GENES = PROJECT / "05_annotation/liftoff_final/cervil.final.liftoff.gff3"
TES = PROJECT / "05_annotation/edta_final/ragtag.scaffold.fasta.mod.EDTA.intact.gff3"

OUT = REPORT / "top_inversions_gene_TE_overlap.tsv"

def chr_to_ragtag(x):
    return x.replace("Chr", "") + "_RagTag" if x.startswith("Chr") else x

def read_gff_intervals(gff, feature_keep=None):
    d = {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip().split("\t")
            if len(p) < 9:
                continue
            chrom, feat = p[0], p[2]
            if feature_keep and feat not in feature_keep:
                continue
            start, end = int(p[3]), int(p[4])
            d.setdefault(chrom, []).append((start, end, feat, p[8]))
    for chrom in d:
        d[chrom].sort()
    return d

def overlap_features(chrom, start, end, intervals):
    hits = []
    for s, e, feat, attr in intervals.get(chrom, []):
        if e < start:
            continue
        if s > end:
            break
        ov = max(0, min(end, e) - max(start, s) + 1)
        if ov > 0:
            hits.append((s, e, feat, ov, attr))
    return hits

genes = read_gff_intervals(GENES, feature_keep={"gene"})
tes = read_gff_intervals(TES, feature_keep=None)

rows = []
with open(EVENTS) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if row["event_type"] != "INV":
            continue
        if row["confidence_class"] != "high_confidence_candidate":
            continue

        qchr = chr_to_ragtag(row["qry_chr"])
        qstart = int(row["qry_start"])
        qend = int(row["qry_end"])
        if qstart > qend:
            qstart, qend = qend, qstart

        gene_hits = overlap_features(qchr, qstart, qend, genes)
        te_hits = overlap_features(qchr, qstart, qend, tes)

        inv_len = int(row["ref_len_bp"])
        te_bp = sum(h[3] for h in te_hits)
        te_frac = te_bp / max(1, qend - qstart + 1)

        rows.append({
            "ref_chr": row["ref_chr"],
            "ref_start": row["ref_start"],
            "ref_end": row["ref_end"],
            "ref_len_bp": row["ref_len_bp"],
            "qry_chr": qchr,
            "qry_start": qstart,
            "qry_end": qend,
            "qry_len_bp": row["qry_len_bp"],
            "n_genes_overlap": len(gene_hits),
            "n_TE_features_overlap": len(te_hits),
            "TE_bp_overlap": te_bp,
            "TE_fraction_query": round(te_frac, 4),
            "near_ragtag_join_10kb": row["near_ragtag_join_10kb"],
            "confidence_class": row["confidence_class"],
        })

rows.sort(key=lambda x: int(x["ref_len_bp"]), reverse=True)

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"Wrote: {OUT}")
