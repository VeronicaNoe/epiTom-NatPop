#!/usr/bin/env bash
set -euo pipefail

PUB="/mnt/disk2/vibanez/otherAnalysis/00_public-data/genotypes/genotypes"
IN_SV="${PUB}/706.sv.vcf.gz"

OUTDIR="/mnt/disk2/vibanez/otherAnalysis/12_getting-graphpanSV-sample-names/11_sv-final"
MAP="/mnt/disk2/vibanez/otherAnalysis/12_getting-graphpanSV-sample-names/05_matches/graphpan_to_my.final_rename_map.filtered.tsv"

mkdir -p "$OUTDIR"

KEEP_OLD="${OUTDIR}/keep.graphpan.samples.txt"
SV_SUB_OLD="${OUTDIR}/graphpanSV.selected.oldNames.vcf.gz"
SV_SUB_NEW="${OUTDIR}/graphpanSV.selected.renamed.vcf.gz"

# ------------------------------------------------------------
# checks
# ------------------------------------------------------------
[[ -f "$IN_SV" ]] || { echo "[ERROR] Missing $IN_SV" >&2; exit 1; }
[[ -f "$MAP"   ]] || { echo "[ERROR] Missing $MAP" >&2; exit 1; }

# old graphpan sample names to keep
cut -f1 "$MAP" > "$KEEP_OLD"

# duplicated final names are not allowed
dup_new=$(cut -f2 "$MAP" | sort | uniq -d || true)
if [[ -n "$dup_new" ]]; then
  echo "[ERROR] Duplicated new sample names in rename map:" >&2
  echo "$dup_new" >&2
  exit 1
fi

# check all old names exist in the SV VCF
bcftools query -l "$IN_SV" | sort > "${OUTDIR}/all_sv_samples.txt"
sort "$KEEP_OLD" > "${OUTDIR}/keep_old.sorted.txt"

missing_old=$(comm -23 "${OUTDIR}/keep_old.sorted.txt" "${OUTDIR}/all_sv_samples.txt" || true)
if [[ -n "$missing_old" ]]; then
  echo "[ERROR] These graphpan samples from the map are not in $IN_SV:" >&2
  echo "$missing_old" >&2
  exit 1
fi

# ------------------------------------------------------------
# 1) subset only the selected graphpan samples
# ------------------------------------------------------------
bcftools view \
  -S "$KEEP_OLD" \
  -Oz -o "$SV_SUB_OLD" \
  "$IN_SV"

bcftools index -t "$SV_SUB_OLD"

# ------------------------------------------------------------
# 2) rename selected samples using the 2-column map
# ------------------------------------------------------------
bcftools reheader \
  -s "$MAP" \
  -o "$SV_SUB_NEW" \
  "$SV_SUB_OLD"

bcftools index -t "$SV_SUB_NEW"

# ------------------------------------------------------------
# 3) quick checks
# ------------------------------------------------------------
echo "[INFO] Final sample count:"
bcftools query -l "$SV_SUB_NEW" | wc -l

echo "[INFO] Final sample names:"
bcftools query -l "$SV_SUB_NEW"

echo "[DONE]"
echo "$SV_SUB_NEW"
