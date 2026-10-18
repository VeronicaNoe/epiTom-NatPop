#!/usr/bin/env python3
import csv
import os

DMR_IN = "/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/af_DMR-epi-rate/03.2_gain-loss_rate_density_per_row.tsv"
SMALL_IN = "/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/aa_results/cervil-SNPs-ONT.GQ20_DP5.classified.tsv"
SV_IN = "/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/aa_results/cervil-SV-ONT.SUP5.classified.tsv"
OUTROOT = "/mnt/disk2/vibanez/otherAnalysis/16_pedigree-reference-base-analysis/ag_ONT-DMR-genetic-vars-comparison"

BEDDIR = os.path.join(OUTROOT, "01_bed")
os.makedirs(BEDDIR, exist_ok=True)

INF_CLASSES = {"fixed_G0", "fixed_G5", "seg_G0", "seg_G5"}

def to_int(x):
    return int(float(x))

# -------------------------
# DMR BED
# BED columns:
# chrom, start0, end, dmr_id, DMR, change_type, dir_type, ancestral_state, derived_state
# -------------------------
with open(DMR_IN) as f, open(os.path.join(BEDDIR, "dmrs.bed"), "w") as out:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        chrom = row["chr"]
        start1 = to_int(row["start"])
        end1 = to_int(row["end"])
        start0 = start1 - 1
        dmr_type = row["DMR"]
        change_type = row["change_type"]
        dir_type = row["dir_type"]
        ancestral = row["ancestral_state"]
        derived = row["derived_state"]

        dmr_id = f"{chrom}:{start1}-{end1}:{dmr_type}:{change_type}:{dir_type}"
        out.write(
            "\t".join([
                chrom, str(start0), str(end1), dmr_id,
                dmr_type, change_type, dir_type, ancestral, derived
            ]) + "\n"
        )

# -------------------------
# Small variants BED
# BED columns:
# chrom, start0, end, var_id, var_subtype, var_class, n_g0, n_g5, ref, alt
# -------------------------
with open(SMALL_IN) as f, \
     open(os.path.join(BEDDIR, "smallvars.all.bed"), "w") as out_all, \
     open(os.path.join(BEDDIR, "smallvars.informative.bed"), "w") as out_inf:

    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        chrom = row["chrom"]
        start0 = to_int(row["start"])
        end = to_int(row["end"])
        pos = to_int(row["pos"])
        ref = row["ref"]
        alt = row["alt"]
        var_subtype = row["variant_type"]
        var_class = row["class"]
        n_g0 = row["n_g0"]
        n_g5 = row["n_g5"]

        var_id = f"{chrom}:{pos}:{ref}>{alt}"

        rec = "\t".join([
            chrom, str(start0), str(end), var_id,
            var_subtype, var_class, n_g0, n_g5, ref, alt
        ]) + "\n"

        out_all.write(rec)
        if var_class in INF_CLASSES:
            out_inf.write(rec)

# -------------------------
# SV BED
# BED columns:
# chrom, start0, end, var_id, var_subtype, var_class, n_g0, n_g5, ref, alt
# ref/alt are '.' here just to keep same column count as smallvars
# -------------------------
with open(SV_IN) as f, \
     open(os.path.join(BEDDIR, "sv.all.bed"), "w") as out_all, \
     open(os.path.join(BEDDIR, "sv.informative.bed"), "w") as out_inf:

    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        chrom = row["chrom"]
        start0 = to_int(row["start"])
        end = to_int(row["end"])
        pos = to_int(row["pos"])
        var_subtype = row["svtype"]
        var_class = row["class"]
        n_g0 = row["n_g0"]
        n_g5 = row["n_g5"]

        var_id = f"{chrom}:{pos}:{var_subtype}:{end}"

        rec = "\t".join([
            chrom, str(start0), str(end), var_id,
            var_subtype, var_class, n_g0, n_g5, ".", "."
        ]) + "\n"

        out_all.write(rec)
        if var_class in INF_CLASSES:
            out_inf.write(rec)

print("BED files written to:", BEDDIR)
