#!/usr/bin/env bash
set -euo pipefail

OUTDIR="/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/ad_heritability-merged"
mkdir -p "$OUTDIR"

AUDIT_TSV="${OUTDIR}/00.0_file-audit_all-markers.tsv"

dmr_classes=("C-DMR" "CG-DMR")
types_h2=("sig" "nonSig" "nonSig-GIF")
types_beta=("sig")
chrs=()
for i in $(seq -w 1 12); do
  chrs+=("ch${i}")
done

SNP_ROOT="/mnt/disk2/vibanez/10_data-analysis/Fig3/aa_GWAS-DMRs/bd_results"
TIP_ROOT="/mnt/disk2/vibanez/otherAnalysis/11_DMR-GWAS-TIPs/ab_results"
SV_ROOT="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/ad_results"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

echo -e "marker_set\tdmr_class\tchr\ttype\texpected_reml\texpected_mqtl\tn_reml\tn_mqtl\tn_common\tn_reml_only\tn_mqtl_only" > "$AUDIT_TSV"

get_files_nested() {
  local root="$1" dmr="$2" chr="$3" type="$4" ext="$5"
  local dir="${root}/${dmr}/${chr}/${type}"
  if [[ -d "$dir" ]]; then
    find "$dir" -maxdepth 1 -type f -name "*.${ext}" -printf "%f\n" | sed "s/\.${ext}$//" | sort -u
  fi
}

get_files_flat() {
  local root="$1" dmr="$2" chr="$3" type="$4" ext="$5"
  local dir="${root}/${type}"
  if [[ -d "$dir" ]]; then
    find "$dir" -maxdepth 1 -type f -name "${chr}_${dmr}_*.${ext}" -printf "%f\n" | sed "s/\.${ext}$//" | sort -u
  fi
}

audit_one() {
  local marker="$1" layout="$2" root="$3" dmr="$4" chr="$5" type="$6" expected_reml="$7" expected_mqtl="$8"

  local reml_ids="${tmpdir}/${marker}_${dmr}_${chr}_${type}.reml.ids"
  local mqtl_ids="${tmpdir}/${marker}_${dmr}_${chr}_${type}.mqtl.ids"

  if [[ "$layout" == "flat" ]]; then
    get_files_flat   "$root" "$dmr" "$chr" "$type" "reml" > "$reml_ids"
    get_files_flat   "$root" "$dmr" "$chr" "$type" "mQTL" > "$mqtl_ids"
  else
    get_files_nested "$root" "$dmr" "$chr" "$type" "reml" > "$reml_ids"
    get_files_nested "$root" "$dmr" "$chr" "$type" "mQTL" > "$mqtl_ids"
  fi

  local n_reml n_mqtl n_common n_reml_only n_mqtl_only
  n_reml=$(wc -l < "$reml_ids")
  n_mqtl=$(wc -l < "$mqtl_ids")
  n_common=$(comm -12 "$reml_ids" "$mqtl_ids" | wc -l)
  n_reml_only=$(comm -23 "$reml_ids" "$mqtl_ids" | wc -l)
  n_mqtl_only=$(comm -13 "$reml_ids" "$mqtl_ids" | wc -l)

  echo -e "${marker}\t${dmr}\t${chr}\t${type}\t${expected_reml}\t${expected_mqtl}\t${n_reml}\t${n_mqtl}\t${n_common}\t${n_reml_only}\t${n_mqtl_only}" >> "$AUDIT_TSV"
}

for dmr in "${dmr_classes[@]}"; do
  for chr in "${chrs[@]}"; do
    for type in "${types_h2[@]}"; do
      audit_one "SNP" "flat"   "$SNP_ROOT" "$dmr" "$chr" "$type" "TRUE"  "FALSE"
      audit_one "TIP" "nested" "$TIP_ROOT" "$dmr" "$chr" "$type" "TRUE"  "FALSE"
      audit_one "SV"  "nested" "$SV_ROOT"  "$dmr" "$chr" "$type" "TRUE"  "FALSE"
    done
    for type in "${types_beta[@]}"; do
      audit_one "SNP" "flat"   "$SNP_ROOT" "$dmr" "$chr" "$type" "FALSE" "TRUE"
      audit_one "TIP" "nested" "$TIP_ROOT" "$dmr" "$chr" "$type" "FALSE" "TRUE"
      audit_one "SV"  "nested" "$SV_ROOT"  "$dmr" "$chr" "$type" "FALSE" "TRUE"
    done
  done
done

echo "Saved audit table: $AUDIT_TSV"
echo

echo "Jobs missing expected .reml:"
awk -F'\t' 'NR==1 || ($5=="TRUE" && $7==0)' "$AUDIT_TSV"

echo
echo "Jobs missing expected .mQTL:"
awk -F'\t' 'NR==1 || ($6=="TRUE" && $8==0)' "$AUDIT_TSV"

echo
echo "Jobs with mismatched IDs between .reml and .mQTL:"
awk -F'\t' 'NR==1 || ($10>0 || $11>0)' "$AUDIT_TSV"

echo
echo "Quick summary:"
awk -F'\t' '
NR==1 {next}
{
  key=$1 FS $2 FS $4
  jobs[key]++
  reml_jobs[key]+=($7>0)
  mqtl_jobs[key]+=($8>0)
  total_reml[key]+=$7
  total_mqtl[key]+=$8
}
END{
  print "marker_set\tdmr_class\ttype\tjobs\tjobs_with_reml\tjobs_with_mqtl\ttotal_reml\ttotal_mqtl"
  for (k in jobs) {
    split(k,a,FS)
    print a[1] FS a[2] FS a[3] FS jobs[k] FS reml_jobs[k] FS mqtl_jobs[k] FS total_reml[k] FS total_mqtl[k]
  }
}' "$AUDIT_TSV" | sort
