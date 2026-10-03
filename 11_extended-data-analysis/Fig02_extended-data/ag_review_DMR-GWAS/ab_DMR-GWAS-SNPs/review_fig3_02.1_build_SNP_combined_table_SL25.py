#!/users/bioinfo/vibanez/anaconda3/bin/python

import os
import re
import gzip
import glob
import argparse
import numpy as np
from scipy.stats import chi2


# =========================
# Arguments
# =========================
parser = argparse.ArgumentParser(
    description="Build SNP combined Fisher table using only ps.gz files matching current sig/*.mQTL DMRs."
)
parser.add_argument("dmr_type", choices=["C-DMR", "CG-DMR"])
args = parser.parse_args()

DMR_TYPE = args.dmr_type


# =========================
# Paths
# =========================
# Folder containing the CURRENT significant DMR list as *.mQTL
SIG_MQTL_DIR = "/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results/sig"

# Folder containing the matching ps.gz files
# In your current setup, the selected ps.gz were extracted into the same sig folder
PS_DIR = "/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results/sig"

# SNP BIM in SL2.5
BIM_FILE = "/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/ba_markers/SNPs.LD.bim"

# SL2.5 FAI
SL25_FAI = "/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/ITAG2.4_genomic.fasta.fai"

# Output folder
RESULTS_DIR = "/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/ab_DMR-GWAS-SNPs"
os.makedirs(RESULTS_DIR, exist_ok=True)

COMBINED_TSV = os.path.join(RESULTS_DIR, f"combined_fisher_SNP_{DMR_TYPE}.tsv")

P_FLOOR = 1e-300
PRIMARY_CHROMS = [f"{i:02d}" for i in range(1, 13)]
BIN_SIZE = 1_000_000


# =========================
# Helpers
# =========================
def normalize_chr(chrn):
    s = str(chrn).strip()
    s = re.sub(r"^SL2\.50ch", "", s, flags=re.IGNORECASE)
    s = re.sub(r"^(chr|ch)", "", s, flags=re.IGNORECASE)
    m = re.search(r"(\d+)$", s)
    if m is None:
        return None
    return f"{int(m.group(1)):02d}"


def load_chr_sizes_from_fai(fai_file):
    chr_sizes = {}

    with open(fai_file, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            chrom = normalize_chr(parts[0])
            if chrom is None:
                continue

            if chrom in PRIMARY_CHROMS:
                chr_sizes[chrom] = int(parts[1])

    missing = [c for c in PRIMARY_CHROMS if c not in chr_sizes]
    if missing:
        raise ValueError(f"Missing chromosomes in FAI: {missing}")

    return chr_sizes


def linear_genome_pos(chrom, pos, cumulative_lengths):
    return cumulative_lengths[chrom] + pos


def safe_pvalue(p):
    if not np.isfinite(p):
        return None
    if p <= 0:
        return P_FLOOR
    if p > 1:
        return None
    return p


def collect_selected_dmr_ids(sig_dir, dmr_type):
    pattern = os.path.join(sig_dir, f"*{dmr_type}*.mQTL")
    files = sorted(glob.glob(pattern))
    dmr_ids = [re.sub(r"\.mQTL$", "", os.path.basename(x)) for x in files]
    return sorted(set(dmr_ids))


def select_ps_files_from_sig_dir(sig_dmr_ids, ps_dir):
    selected = []
    missing = []

    for dmr_id in sig_dmr_ids:
        path = os.path.join(ps_dir, f"{dmr_id}.ps.gz")
        if os.path.exists(path):
            selected.append(path)
        else:
            missing.append(dmr_id)

    return selected, missing


def load_snp_pos_map_from_bim(bim_file):
    snp_pos_map = {}
    snp_ids = []

    with open(bim_file, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue

            chrom = normalize_chr(parts[0])
            snp_id = parts[1]

            try:
                pos = int(parts[3])
            except ValueError:
                continue

            if chrom is None or chrom not in PRIMARY_CHROMS:
                continue

            snp_pos_map[snp_id] = (chrom, pos)
            snp_ids.append(snp_id)

    return snp_ids, snp_pos_map


def save_combined_results(
    out_tsv,
    snp_ids,
    sum_log_p,
    file_count,
    snp_pos_map,
    cumulative_lengths,
    bin_size
):
    valid = file_count > 0
    combined_chisq = np.zeros(len(snp_ids), dtype=np.float64)
    combined_pvalue = np.ones(len(snp_ids), dtype=np.float64)

    combined_chisq[valid] = -2.0 * sum_log_p[valid]
    combined_pvalue[valid] = chi2.sf(combined_chisq[valid], 2 * file_count[valid])

    n_written = 0
    with open(out_tsv, "w") as out:
        out.write("key\tchrom\tpos\tlinear_pos\tx_bin\tcombined_chisq\tcombined_pvalue\tfile_count\n")

        for i, snp_id in enumerate(snp_ids):
            if not valid[i]:
                continue
            if snp_id not in snp_pos_map:
                continue

            chrom, pos = snp_pos_map[snp_id]
            linpos = linear_genome_pos(chrom, pos, cumulative_lengths)
            x_bin = linpos / bin_size

            out.write(
                f"{snp_id}\t{chrom}\t{pos}\t{linpos}\t{x_bin:.6f}\t"
                f"{combined_chisq[i]:.6f}\t{combined_pvalue[i]:.6e}\t{file_count[i]}\n"
            )
            n_written += 1

    print(f"Wrote {n_written} rows to: {out_tsv}")


# =========================
# Genome layout
# =========================
chr_sizes = load_chr_sizes_from_fai(SL25_FAI)
chrom_order = [f"{i:02d}" for i in range(1, 13)]

cumulative_lengths = {}
cum_len = 0
for chrom in chrom_order:
    cumulative_lengths[chrom] = cum_len
    cum_len += chr_sizes[chrom]

print("Chromosome sizes loaded from SL2.5 FAI:")
for chrom in chrom_order:
    print(f"{chrom}\t{chr_sizes[chrom]}")


# =========================
# 1) DMR IDs from current sig/*.mQTL
# =========================
sig_dmr_ids = collect_selected_dmr_ids(SIG_MQTL_DIR, DMR_TYPE)
print(f"Selected DMRs from sig/*.mQTL: {len(sig_dmr_ids)}")

if len(sig_dmr_ids) == 0:
    raise ValueError(f"No sig/*.mQTL files found for {DMR_TYPE} in {SIG_MQTL_DIR}")


# =========================
# 2) Match to ps.gz
# =========================
input_files, missing_ps = select_ps_files_from_sig_dir(sig_dmr_ids, PS_DIR)

print(f"Matched ps.gz files: {len(input_files)}")
print(f"Missing ps.gz for selected DMRs: {len(missing_ps)}")
if len(missing_ps) > 0:
    print("First 10 missing:")
    for x in missing_ps[:10]:
        print(f"  {x}")

if len(input_files) == 0:
    raise ValueError("No matching ps.gz files found for the selected sig DMRs")


# =========================
# 3) SNP positions from BIM
# =========================
snp_ids, snp_pos_map = load_snp_pos_map_from_bim(BIM_FILE)
print(f"Loaded SNP positions from BIM: {len(snp_ids)}")

snp_index = {snp_id: i for i, snp_id in enumerate(snp_ids)}
sum_log_p = np.zeros(len(snp_ids), dtype=np.float64)
file_count = np.zeros(len(snp_ids), dtype=np.int32)


# =========================
# 4) Combine Fisher across selected ps.gz
# =========================
for file_idx, file_path in enumerate(input_files, start=1):
    if file_idx % 100 == 0 or file_idx == 1 or file_idx == len(input_files):
        print(f"Processing file {file_idx}/{len(input_files)}: {os.path.basename(file_path)}")

    with gzip.open(file_path, "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue

            snp_id = parts[0]
            idx = snp_index.get(snp_id)
            if idx is None:
                continue

            try:
                raw_p = float(parts[3]) if parts[3] != "nan" else np.nan
            except ValueError:
                continue

            p = safe_pvalue(raw_p)
            if p is None:
                continue

            sum_log_p[idx] += np.log(p)
            file_count[idx] += 1


# =========================
# 5) Save output
# =========================
save_combined_results(
    COMBINED_TSV,
    snp_ids,
    sum_log_p,
    file_count,
    snp_pos_map,
    cumulative_lengths,
    BIN_SIZE
)

print(f"Done: {COMBINED_TSV}")
