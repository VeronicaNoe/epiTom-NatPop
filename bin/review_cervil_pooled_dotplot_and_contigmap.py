#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# =========================
# Paths
# =========================
PROJECT = Path("/work/project/ijpb/vibanez/14_cervil-project/07_cervil-denovo-assembly")
REPORT  = PROJECT / "06_report"

REF_FAI = Path("/home/vibanez/work/01_mapping/01.0_ref/SL5.0.fasta.fai")

DATASETS = {
    "hifiasm_raw": {
        "paf": REPORT / "hifiasm_raw.vs_SL5.paf",
        "query_fai": PROJECT / "01_hifiasm/cervil.primary.fasta.fai",
        "plot_query_fai": PROJECT / "01_hifiasm/cervil.primary.fasta.fai",
        "force_plot_query_order": False,
    },
    "dorado_polish": {
        "paf": REPORT / "dorado_polish.vs_SL5.paf",
        "query_fai": PROJECT / "02_polish/cervil.polished.fasta.fai",
        "plot_query_fai": PROJECT / "02_polish/cervil.polished.fasta.fai",
        "force_plot_query_order": False,
    },
    "ragtag": {
        "paf": REPORT / "ragtag.vs_SL5.paf",
        "query_fai": PROJECT / "03_ragtag/ragtag.scaffold.fasta.fai",

        # Plot only the 12 real chromosome-like RagTag scaffolds.
        # Tables still use the full ragtag.scaffold.fasta.fai.
        "plot_query_fai": REPORT / "cervil.ragtag.chr_like.fasta.fai",
        "force_plot_query_order": True,
    },
}

OUTDIR = REPORT / "pooled_qc"
DOTDIR = OUTDIR / "dotplots"
TABDIR = OUTDIR / "tables"
DOTDIR.mkdir(parents=True, exist_ok=True)
TABDIR.mkdir(parents=True, exist_ok=True)

# =========================
# Plot filter: relaxed enough to show genome structure
# =========================
PLOT_FILTER = {
    "min_aln": 1000,
    "min_mapq": 0,
    "min_identity": 0,
    "drop_secondary": False,
}

# =========================
# Table filters
# =========================
TABLE_FILTERS = {
    "strict": {
        "min_aln": 10000,
        "min_mapq": 20,
        "min_identity": 0.80,
        "drop_secondary": True,
        "low_aligned_frac": 0.20,
        "unique_main_frac": 0.80,
    },
    "relaxed": {
        "min_aln": 1000,
        "min_mapq": 0,
        "min_identity": 0.75,
        "drop_secondary": False,
        "low_aligned_frac": 0.10,
        "unique_main_frac": 0.80,
    },
}

# Only real SL5 chromosomes for plotting.
SL5_CHR_ORDER = [str(i) for i in range(1, 13)]

# Only real RagTag chromosome-like scaffolds for the RagTag plot.
RAGTAG_CHR_ORDER = [f"{i}_RagTag" for i in range(1, 13)]


# =========================
# Helpers
# =========================
def read_fai(fai):
    rows = []
    with open(fai) as f:
        for i, line in enumerate(f):
            p = line.rstrip().split("\t")
            rows.append({"name": p[0], "length": int(p[1]), "order": i})
    df = pd.DataFrame(rows)
    df["offset"] = df["length"].cumsum() - df["length"]
    return df


def subset_fai_by_names(fai_df, names):
    names = [x for x in names if x in set(fai_df["name"])]
    out = fai_df[fai_df["name"].isin(names)].copy()
    out["name"] = pd.Categorical(out["name"], categories=names, ordered=True)
    out = out.sort_values("name").copy()
    out["name"] = out["name"].astype(str)
    out["order"] = range(len(out))
    out["offset"] = out["length"].cumsum() - out["length"]
    return out


def read_paf(paf):
    rows = []

    with open(paf) as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 12:
                continue

            tp = "NA"
            for tag in p[12:]:
                if tag.startswith("tp:A:"):
                    tp = tag.split(":")[-1]

            qname = p[0]
            qlen = int(p[1])
            qstart = int(p[2])
            qend = int(p[3])
            strand = p[4]
            tname = p[5]
            tlen = int(p[6])
            tstart = int(p[7])
            tend = int(p[8])
            nmatch = int(p[9])
            alen = int(p[10])
            mapq = int(p[11])
            identity = nmatch / alen if alen > 0 else 0

            rows.append({
                "qname": qname,
                "qlen": qlen,
                "qstart": qstart,
                "qend": qend,
                "strand": strand,
                "tname": tname,
                "tlen": tlen,
                "tstart": tstart,
                "tend": tend,
                "nmatch": nmatch,
                "alen": alen,
                "mapq": mapq,
                "identity": identity,
                "tp": tp,
            })

    return pd.DataFrame(rows)


def apply_filter(df, filt):
    if df.empty:
        return df.copy()

    x = df[
        (df["alen"] >= filt["min_aln"]) &
        (df["mapq"] >= filt["min_mapq"]) &
        (df["identity"] >= filt["min_identity"])
    ].copy()

    if filt["drop_secondary"]:
        x = x[x["tp"] != "S"].copy()

    return x


def build_offsets(names, lengths):
    offsets = {}
    cur = 0
    for name in names:
        offsets[name] = cur
        cur += int(lengths[name])
    return offsets, cur


def union_length(intervals):
    if not intervals:
        return 0

    intervals = sorted(intervals)
    total = 0
    cs, ce = intervals[0]

    for s, e in intervals[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            total += ce - cs
            cs, ce = s, e

    total += ce - cs
    return total


def order_queries_by_main_SL5(df, qfai, tfai):
    allowed_q = set(qfai["name"])
    allowed_t = set(tfai["name"])
    d = df[df["qname"].isin(allowed_q) & df["tname"].isin(allowed_t)].copy()

    qlen = dict(zip(qfai["name"], qfai["length"]))
    target_order = dict(zip(tfai["name"], tfai["order"]))
    target_offset = dict(zip(tfai["name"], tfai["offset"]))

    if d.empty:
        qnames = list(qfai.sort_values("length", ascending=False)["name"])
        odf = pd.DataFrame({
            "qname": qnames,
            "main_SL5_chr": "NA",
            "target_mid_cum": pd.NA,
            "query_len": [qlen[x] for x in qnames],
        })
        return qnames, odf

    d["target_mid_cum"] = d["tname"].map(target_offset) + (d["tstart"] + d["tend"]) / 2

    rows = []
    for qname, g in d.groupby("qname"):
        by_chr = g.groupby("tname")["alen"].sum().sort_values(ascending=False)
        main_chr = by_chr.index[0]
        main = g[g["tname"] == main_chr]

        rows.append({
            "qname": qname,
            "main_SL5_chr": main_chr,
            "main_chr_order": target_order.get(main_chr, 10**9),
            "target_mid_cum": main["target_mid_cum"].median(),
            "query_len": qlen.get(qname, 0),
        })

    odf = pd.DataFrame(rows)
    odf = odf.sort_values(
        ["main_chr_order", "target_mid_cum", "query_len"],
        ascending=[True, True, False]
    ).copy()

    ordered = list(odf["qname"])
    missing = [x for x in qfai["name"] if x not in set(ordered)]
    missing = sorted(missing, key=lambda z: qlen[z], reverse=True)

    if missing:
        miss = pd.DataFrame({
            "qname": missing,
            "main_SL5_chr": "NA",
            "main_chr_order": 10**9,
            "target_mid_cum": pd.NA,
            "query_len": [qlen[x] for x in missing],
        })
        odf = pd.concat([odf, miss], ignore_index=True)

    return ordered + missing, odf


def make_forced_order_table(qnames):
    rows = []
    for q in qnames:
        main_chr = q.replace("_RagTag", "")
        rows.append({
            "qname": q,
            "main_SL5_chr": main_chr,
            "main_chr_order": int(main_chr) if main_chr.isdigit() else 10**9,
            "target_mid_cum": pd.NA,
            "query_len": pd.NA,
        })
    return pd.DataFrame(rows)


def make_contig_table(state, df, qfai, filt_name, filt):
    qlen = dict(zip(qfai["name"], qfai["length"]))
    rows = []

    for qname, length in qlen.items():
        g = df[df["qname"] == qname]

        if g.empty:
            rows.append({
                "state": state,
                "filter": filt_name,
                "query_contig": qname,
                "query_len": length,
                "main_SL5_chr": "NA",
                "aligned_bp_to_main_chr": 0,
                "second_SL5_chr": "NA",
                "aligned_bp_to_second_chr": 0,
                "total_aligned_bp": 0,
                "query_union_aligned_bp": 0,
                "pct_query_aligned": 0,
                "pct_aligned_to_main_chr": 0,
                "n_SL5_chr_hit": 0,
                "status": "unplaced",
            })
            continue

        by_chr = g.groupby("tname")["alen"].sum().sort_values(ascending=False)

        main_chr = by_chr.index[0]
        main_bp = int(by_chr.iloc[0])

        second_chr = by_chr.index[1] if len(by_chr) > 1 else "NA"
        second_bp = int(by_chr.iloc[1]) if len(by_chr) > 1 else 0

        total_bp = int(by_chr.sum())
        intervals = [(int(r.qstart), int(r.qend)) for r in g.itertuples()]
        query_union_bp = union_length(intervals)

        pct_query_aligned = query_union_bp / length if length else 0
        pct_main = main_bp / total_bp if total_bp else 0
        n_chr_hit = int((by_chr >= filt["min_aln"]).sum())

        if pct_query_aligned < filt["low_aligned_frac"]:
            status = "low_aligned"
        elif pct_main >= filt["unique_main_frac"]:
            status = "unique_chr"
        else:
            status = "multi_chr"

        rows.append({
            "state": state,
            "filter": filt_name,
            "query_contig": qname,
            "query_len": length,
            "main_SL5_chr": main_chr,
            "aligned_bp_to_main_chr": main_bp,
            "second_SL5_chr": second_chr,
            "aligned_bp_to_second_chr": second_bp,
            "total_aligned_bp": total_bp,
            "query_union_aligned_bp": query_union_bp,
            "pct_query_aligned": pct_query_aligned,
            "pct_aligned_to_main_chr": pct_main,
            "n_SL5_chr_hit": n_chr_hit,
            "status": status,
        })

    out = pd.DataFrame(rows)

    status_rank = {
        "unique_chr": 0,
        "multi_chr": 1,
        "low_aligned": 2,
        "unplaced": 3,
    }

    out["status_rank"] = out["status"].map(status_rank).fillna(9)
    out = out.sort_values(
        ["status_rank", "query_len"],
        ascending=[True, False]
    ).drop(columns=["status_rank"])

    return out


def get_plot_query_fai(state, cfg, full_qfai):
    plot_fai = Path(cfg["plot_query_fai"])

    if state == "ragtag":
        if not plot_fai.exists():
            raise FileNotFoundError(plot_fai)

        qfai = read_fai(plot_fai)
        qfai = subset_fai_by_names(qfai, RAGTAG_CHR_ORDER)

        missing = [x for x in RAGTAG_CHR_ORDER if x not in set(qfai["name"])]
        if missing:
            raise ValueError(f"Missing RagTag chr-like scaffolds in {plot_fai}: {missing}")

        return qfai

    if plot_fai.exists():
        return read_fai(plot_fai)

    return full_qfai.copy()


def clean_label(name):
    return str(name).replace("_RagTag", "")


def plot_dotplot(state, df, qfai_plot, tfai_plot, force_query_order=False):
    allowed_q = set(qfai_plot["name"])
    allowed_t = set(tfai_plot["name"])

    d = df[df["qname"].isin(allowed_q) & df["tname"].isin(allowed_t)].copy()

    qlen = dict(zip(qfai_plot["name"], qfai_plot["length"]))
    tlen = dict(zip(tfai_plot["name"], tfai_plot["length"]))

    if force_query_order:
        qnames = list(qfai_plot["name"])
        order_table = make_forced_order_table(qnames)
    else:
        qnames, order_table = order_queries_by_main_SL5(d, qfai_plot, tfai_plot)

    tnames = list(tfai_plot["name"])

    qoff, qtotal = build_offsets(qnames, qlen)
    toff, ttotal = build_offsets(tnames, tlen)

    plus_segments = []
    minus_segments = []

    for r in d.itertuples():
        if r.tname not in toff or r.qname not in qoff:
            continue

        x1 = toff[r.tname] + int(r.tstart)
        x2 = toff[r.tname] + int(r.tend)

        if r.strand == "+":
            y1 = qoff[r.qname] + int(r.qstart)
            y2 = qoff[r.qname] + int(r.qend)
            plus_segments.append([(x1, y1), (x2, y2)])
        else:
            y1 = qoff[r.qname] + int(r.qend)
            y2 = qoff[r.qname] + int(r.qstart)
            minus_segments.append([(x1, y1), (x2, y2)])

    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)

    if plus_segments:
        ax.add_collection(LineCollection(
            plus_segments,
            colors="black",
            linewidths=0.30,
            alpha=0.70
        ))

    if minus_segments:
        ax.add_collection(LineCollection(
            minus_segments,
            colors="red",
            linewidths=0.30,
            alpha=0.70
        ))

    # Chromosome/contig boundaries
    for x in [toff[t] for t in tnames]:
        ax.axvline(x, color="grey", linewidth=0.25, alpha=0.35)

    if force_query_order or len(qnames) <= 30:
        for q in qnames:
            ax.axhline(qoff[q], color="grey", linewidth=0.25, alpha=0.35)
    else:
        # For contig-level assemblies, draw only broad group boundaries.
        seen = set()
        for _, r in order_table.iterrows():
            main_chr = r["main_SL5_chr"]
            q = r["qname"]
            if main_chr in seen or main_chr == "NA":
                continue
            if q not in qoff:
                continue
            seen.add(main_chr)
            ax.axhline(qoff[q], color="grey", linewidth=0.25, alpha=0.35)

    ax.set_xlim(0, ttotal)
    ax.set_ylim(0, qtotal)

    ax.set_xlabel("SL5 cumulative coordinate")
    if state == "ragtag":
        ax.set_ylabel("Cervil RagTag chromosome-like scaffolds")
    else:
        ax.set_ylabel("Cervil assembly cumulative coordinate\n(contigs ordered by main SL5 hit)")

    ax.set_title(f"{state} vs SL5")

    # Clean chr labels: only 1-12
    xticks = [toff[t] + tlen[t] / 2 for t in tnames]
    xlabels = [clean_label(t) for t in tnames]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=90, fontsize=7)

    if force_query_order or len(qnames) <= 30:
        yticks = [qoff[q] + qlen[q] / 2 for q in qnames]
        ylabels = [clean_label(q) for q in qnames]
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, fontsize=7)

    pdf = DOTDIR / f"{state}.vs_SL5.cumulative.clear.pdf"
    svg = DOTDIR / f"{state}.vs_SL5.cumulative.clear.svg"
    order_out = TABDIR / f"{state}.plot_query_order.tsv"

    fig.savefig(pdf)
    fig.savefig(svg)
    plt.close(fig)

    order_table.to_csv(order_out, sep="\t", index=False)

    print(f"[dotplot] {pdf}")
    print(f"[dotplot] {svg}")
    print(f"[plot order] {order_out}")
    print(f"[plot alignments] {state}: {len(d)}")


def main():
    tfai_all = read_fai(REF_FAI)

    # Plot only real SL5 chromosomes 1-12.
    tfai_plot = subset_fai_by_names(tfai_all, SL5_CHR_ORDER)

    missing_ref = [x for x in SL5_CHR_ORDER if x not in set(tfai_plot["name"])]
    if missing_ref:
        raise ValueError(f"Missing SL5 chromosomes in {REF_FAI}: {missing_ref}")

    all_tables = []
    filter_summaries = []

    for state, cfg in DATASETS.items():
        paf = cfg["paf"]
        qfai_path = cfg["query_fai"]

        if not paf.exists():
            raise FileNotFoundError(paf)
        if not qfai_path.exists():
            raise FileNotFoundError(qfai_path)

        print(f"\n[state] {state}")
        print(f"[paf] {paf}")
        print(f"[query fai] {qfai_path}")

        qfai_full = read_fai(qfai_path)
        qfai_plot = get_plot_query_fai(state, cfg, qfai_full)

        raw = read_paf(paf)
        print(f"[raw PAF hits] {len(raw)}")

        # Plot uses relaxed visual filter.
        plot_df = apply_filter(raw, PLOT_FILTER)
        print(f"[plot PAF hits kept before chr subset] {len(plot_df)}")

        filter_summaries.append({
            "state": state,
            "filter": "plot",
            "min_aln": PLOT_FILTER["min_aln"],
            "min_mapq": PLOT_FILTER["min_mapq"],
            "min_identity": PLOT_FILTER["min_identity"],
            "drop_secondary": PLOT_FILTER["drop_secondary"],
            "n_paf_hits": len(plot_df),
            "total_aligned_bp": int(plot_df["alen"].sum()) if not plot_df.empty else 0,
            "n_query_with_hit": int(plot_df["qname"].nunique()) if not plot_df.empty else 0,
            "n_total_query": len(qfai_full),
        })

        plot_dotplot(
            state,
            plot_df,
            qfai_plot,
            tfai_plot,
            force_query_order=cfg.get("force_plot_query_order", False),
        )

        # Tables: strict + relaxed, using full query assembly and all SL5 targets.
        for filt_name, filt in TABLE_FILTERS.items():
            table_df = apply_filter(raw, filt)
            print(f"[{filt_name} PAF hits kept] {len(table_df)}")

            tab = make_contig_table(state, table_df, qfai_full, filt_name, filt)
            outtab = TABDIR / f"{state}.contig_to_SL5.{filt_name}.tsv"
            tab.to_csv(outtab, sep="\t", index=False)
            print(f"[table] {outtab}")

            all_tables.append(tab)

            filter_summaries.append({
                "state": state,
                "filter": filt_name,
                "min_aln": filt["min_aln"],
                "min_mapq": filt["min_mapq"],
                "min_identity": filt["min_identity"],
                "drop_secondary": filt["drop_secondary"],
                "n_paf_hits": len(table_df),
                "total_aligned_bp": int(table_df["alen"].sum()) if not table_df.empty else 0,
                "n_query_with_hit": int(table_df["qname"].nunique()) if not table_df.empty else 0,
                "n_total_query": len(qfai_full),
            })

    all_tab = pd.concat(all_tables, ignore_index=True)
    all_tab.to_csv(
        TABDIR / "all_states.contig_to_SL5.strict_and_relaxed.tsv",
        sep="\t",
        index=False,
    )

    summary = (
        all_tab
        .groupby(["state", "filter", "status"], as_index=False)
        .agg(
            n_contigs=("query_contig", "count"),
            total_bp=("query_len", "sum"),
            median_pct_query_aligned=("pct_query_aligned", "median"),
            median_pct_main_chr=("pct_aligned_to_main_chr", "median"),
        )
    )
    summary.to_csv(
        TABDIR / "contig_to_SL5.status_summary.strict_and_relaxed.tsv",
        sep="\t",
        index=False,
    )

    filt_summary = pd.DataFrame(filter_summaries)
    filt_summary.to_csv(
        TABDIR / "alignment_filter_summary.tsv",
        sep="\t",
        index=False,
    )

    print("\n[done]")
    print(f"Dotplots: {DOTDIR}")
    print(f"Tables:   {TABDIR}")


if __name__ == "__main__":
    main()
