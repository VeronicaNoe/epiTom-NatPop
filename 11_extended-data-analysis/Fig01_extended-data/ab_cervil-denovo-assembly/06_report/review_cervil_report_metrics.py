#!/usr/bin/env python3
import csv
import os
import re
from pathlib import Path

PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")
QCDIR = PROJECT / "04_qc"
OUT = PROJECT / "06_report/assembly_summary.tsv"

ASSEMBLIES = {
    "hifiasm_raw": PROJECT / "01_hifiasm/cervil.primary.fasta",
    "dorado_polish": PROJECT / "02_polish/cervil.polished.fasta",
    "ragtag": PROJECT / "03_ragtag/ragtag.scaffold.fasta",
}

def fasta_stats(fp):
    lengths = []
    gc = 0
    total = 0
    seqlen = 0
    with open(fp) as f:
        for line in f:
            if line.startswith(">"):
                if seqlen > 0:
                    lengths.append(seqlen)
                seqlen = 0
            else:
                s = line.strip().upper()
                seqlen += len(s)
                total += len(s)
                gc += s.count("G") + s.count("C")
        if seqlen > 0:
            lengths.append(seqlen)

    lengths.sort(reverse=True)
    csum = 0
    n50 = l50 = n90 = l90 = 0
    for i, L in enumerate(lengths, start=1):
        csum += L
        if n50 == 0 and csum >= total * 0.5:
            n50, l50 = L, i
        if n90 == 0 and csum >= total * 0.9:
            n90, l90 = L, i

    return {
        "n_seqs": len(lengths),
        "total_bp": total,
        "largest_bp": lengths[0] if lengths else 0,
        "gc_pct": round(100 * gc / total, 3) if total else 0,
        "n50_bp": n50,
        "l50": l50,
        "n90_bp": n90,
        "l90": l90,
    }

def read_quast_report_tsv(fp):
    vals = {}
    with open(fp) as f:
        reader = csv.reader(f, delimiter="\t")
        rows = list(reader)

    if len(rows) < 2:
        return vals

    header = rows[0]
    value_col = 1
    for row in rows[1:]:
        if len(row) > value_col:
            vals[row[0]] = row[value_col]
    return vals

def read_busco_summary(fp):
    txt = Path(fp).read_text()
    out = {
        "busco_complete_pct": "",
        "busco_single_pct": "",
        "busco_duplicated_pct": "",
        "busco_fragmented_pct": "",
        "busco_missing_pct": "",
        "busco_n": "",
    }
    m = re.search(
        r"C:([\d\.]+)%\[S:([\d\.]+)%,D:([\d\.]+)%\],F:([\d\.]+)%,M:([\d\.]+)%,n:(\d+)",
        txt
    )
    if m:
        out["busco_complete_pct"] = m.group(1)
        out["busco_single_pct"] = m.group(2)
        out["busco_duplicated_pct"] = m.group(3)
        out["busco_fragmented_pct"] = m.group(4)
        out["busco_missing_pct"] = m.group(5)
        out["busco_n"] = m.group(6)
    return out

rows = []
for state, fasta in ASSEMBLIES.items():
    row = {"assembly": state}
    row.update(fasta_stats(fasta))

    quast_tsv = QCDIR / f"quast_{state}/report.tsv"
    if quast_tsv.exists():
        q = read_quast_report_tsv(quast_tsv)
        row["quast_genome_fraction_pct"] = q.get("Genome fraction (%)", "")
        row["quast_dup_ratio"] = q.get("Duplication ratio", "")
        row["quast_misassemblies"] = q.get("# misassemblies", "")
        row["quast_Ns_per_100kbp"] = q.get("# N's per 100 kbp", "")
        row["quast_mismatches_per_100kbp"] = q.get("# mismatches per 100 kbp", "")
        row["quast_indels_per_100kbp"] = q.get("# indels per 100 kbp", "")
    else:
        row["quast_genome_fraction_pct"] = ""
        row["quast_dup_ratio"] = ""
        row["quast_misassemblies"] = ""
        row["quast_Ns_per_100kbp"] = ""
        row["quast_mismatches_per_100kbp"] = ""
        row["quast_indels_per_100kbp"] = ""

    busco_txt = QCDIR / f"busco_{state}/short_summary.specific.embryophyta_odb12.busco_{state}.txt"
    if busco_txt.exists():
        row.update(read_busco_summary(busco_txt))
    else:
        row.update({
            "busco_complete_pct": "",
            "busco_single_pct": "",
            "busco_duplicated_pct": "",
            "busco_fragmented_pct": "",
            "busco_missing_pct": "",
            "busco_n": "",
        })

    rows.append(row)

fields = [
    "assembly",
    "n_seqs", "total_bp", "largest_bp", "gc_pct",
    "n50_bp", "l50", "n90_bp", "l90",
    "quast_genome_fraction_pct", "quast_dup_ratio", "quast_misassemblies",
    "quast_Ns_per_100kbp", "quast_mismatches_per_100kbp", "quast_indels_per_100kbp",
    "busco_complete_pct", "busco_single_pct", "busco_duplicated_pct",
    "busco_fragmented_pct", "busco_missing_pct", "busco_n"
]

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow(r)

print(f"Wrote: {OUT}")
