#!/usr/bin/env bash
set -euo pipefail

TAG="$(cat "$1")"
WDIR="/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS"
TDIR="$WDIR/aa_get-target-dmrs"

TARGETS="$TDIR/${TAG}.targets.tsv"
TFAM="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/ba_markers/DMR_general_leaf_LD.tfam"

SV_VCF="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/aa_markers/graphpanSV.sharedSamples.vcf.gz"
SV_COORD="/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/aa_DMR-GWAS-SV_SL25/sv_id2coord_sl25_primary.tsv"

TOP_VAR_SCAN="${TOP_VAR_SCAN:-2}"
MIN_NONMISS=0
MIN_MAF=0
MIN_MAC=1
KEEP_FAIL_COLS="${KEEP_FAIL_COLS:-0}"

EXPAND_WINS=(50000 100000 250000 500000 1000000 2000000 5000000 10000000)

OUTDIR="$WDIR/ab_covariables/SV"
mkdir -p "$OUTDIR"

CANDS="$OUTDIR/${TAG}.cands.tsv"
COVAR_OUT="$OUTDIR/${TAG}__SV.covar"
INCLUDED="$OUTDIR/${TAG}__SV.included_variants.tsv"

ID3="$OUTDIR/${TAG}.tfam_fid_iid_intercept.tsv"
SAMPLES_FILE="$OUTDIR/${TAG}.samples.tfam_order.txt"
MM_COLMAP="$OUTDIR/${TAG}.SV.colmap.tsv"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

ALL_CAND_IDS="$TMPDIR/all_candidate_ids.txt"
ALL_CAND_TSV="$TMPDIR/all_candidates.tsv"
SV_LINES_TSV="$TMPDIR/sv_lines.tsv"

[[ -s "$TARGETS" ]] || { echo "ERROR: missing TARGETS: $TARGETS" >&2; exit 1; }
[[ -s "$TFAM" ]] || { echo "ERROR: missing TFAM: $TFAM" >&2; exit 1; }
[[ -s "$SV_VCF" ]] || { echo "ERROR: missing SV_VCF: $SV_VCF" >&2; exit 1; }
[[ -s "${SV_VCF}.tbi" ]] || { echo "ERROR: missing VCF index: ${SV_VCF}.tbi" >&2; exit 1; }
[[ -s "$SV_COORD" ]] || { echo "ERROR: missing SV_COORD: $SV_COORD" >&2; exit 1; }

awk 'BEGIN{OFS="\t"}{print $1,$2,1}' "$TFAM" > "$ID3"
awk '{print $2}' "$TFAM" > "$SAMPLES_FILE"

make_colmap() {
  local vcf="$1" out="$2"
  local header
  header="$(bcftools view -h "$vcf" | awk '/^#CHROM/{print; exit}')"
  awk -v hdr="$header" '
    BEGIN{
      FS=OFS="\t";
      n=split(hdr,H,"\t");
      for(i=1;i<=n;i++) col[H[i]]=i;
    }
    function strip_suffix(s){ sub(/_[0-9]+$/, "", s); return s; }
    {
      s=$1;
      if(s in col){ print s,col[s]; next; }
      s2=strip_suffix(s);
      if(s2 in col){ print $1,col[s2]; next; }
      print $1,0;
    }' "$SAMPLES_FILE" > "$out"
}
make_colmap "$SV_VCF" "$MM_COLMAP"

norm_sl25_chr() {
  local raw="$1"
  raw="${raw//$'\r'/}"
  case "$raw" in
    SL2.50ch*) echo "$raw" ;;
    SL2.5ch*)  echo "${raw/SL2.5ch/SL2.50ch}" ;;
    ch[0-9][0-9]) printf "SL2.50%s\n" "$raw" ;;
    chr[0-9]*)
      local n="${raw#chr}"
      printf "SL2.50ch%02d\n" "$((10#$n))"
      ;;
    [0-9]*)
      printf "SL2.50ch%02d\n" "$((10#$raw))"
      ;;
    *)
      echo "$raw"
      ;;
  esac
}

variant_to_dosage_col() {
  local variant_line="$1"
  local colmap="$2"

  printf '%s\n' "$variant_line" | awk -v colmap="$colmap" '
    BEGIN{
      FS=OFS="\t";
      while((getline<colmap)>0){
        samp[++n]=$1;
        idx[$1]=$2;
      }
      close(colmap);
    }
    {
      for(i=1;i<=n;i++){
        c=idx[samp[i]];

        if(c==0){
          print -9;
          continue;
        }

        gtfield=$c;
        split(gtfield,a,":");
        g=a[1];

        if(g=="" || g ~ /\./){
          print -9;
          continue;
        }

        m=split(g,al,/[\/|]/);
        dose=0;
        bad=0;
        for(j=1;j<=m;j++){
          if(al[j] !~ /^[0-9]+$/){
            bad=1;
            break;
          }
          if((al[j]+0)!=0) dose+=1;
        }

        if(bad) print -9;
        else    print dose;
      }
    }'
}

qc_dose() {
  local dose_file="$1"
  awk -v minN="$MIN_NONMISS" -v minMAF="$MIN_MAF" -v minMAC="$MIN_MAC" '
    $1!=-9 {n++; alt+=$1}
    END{
      if(n==0){print "FAIL 0 NA 0"; exit}
      totA=2*n;
      ref=totA-alt;
      mac = (alt<ref ? alt : ref);
      maf = mac/totA;
      ok=1;
      if(n < minN) ok=0;
      if(maf < minMAF) ok=0;
      if(minMAC>0 && mac < minMAC) ok=0;
      if(ok) printf("OK %d %.10g %d\n", n, maf, mac);
      else   printf("FAIL %d %.10g %d\n", n, maf, mac);
    }' "$dose_file"
}

pick_candidates_sv() {
  local coord_tsv="$1" dmr_chr="$2" dmr_pos="$3" k="$4"
  local w tmp n

  for w in "${EXPAND_WINS[@]}"; do
    tmp="$(mktemp)"
    awk -v FS="\t" -v OFS="\t" -v chr="$dmr_chr" -v pos="$dmr_pos" -v w="$w" '
      $4==chr {
        sv_id=$1; sv_chr=$4; sv_start=$5+0; sv_end=$6+0;
        if(pos < sv_start)      dist = sv_start - pos;
        else if(pos > sv_end)   dist = pos - sv_end;
        else                    dist = 0;
        if(dist <= w) print sv_id, sv_chr, sv_start, sv_end, dist;
      }' "$coord_tsv" \
      | sort -t$'\t' -k5,5n -k3,3n -k4,4n \
      | head -n "$k" \
      | awk -v OFS="\t" -v w="$w" '{print $1,$2,$3,$4,$5,w}' > "$tmp"

    n=$(wc -l < "$tmp" | tr -d ' ')
    if [[ "$n" -gt 0 ]]; then
      cat "$tmp"
      rm -f "$tmp"
      return 0
    fi
    rm -f "$tmp"
  done
  return 1
}

get_cached_vline() {
  local id="$1"
  local cache="$2"
  awk -F'\t' -v OFS='\t' -v id="$id" '
    $1==id{
      for(i=2;i<=NF;i++){
        printf "%s", $i
        if(i<NF) printf OFS
      }
      printf "\n"
      exit
    }' "$cache"
}

echo -e "tag\tcov_mm\tqtl_rank\tdmr_id\tdmr_chr\tdmr_pos\tcand_rank\tsv_id\tsv_chr\tsv_start\tsv_end\tdist_bp\tused_win_bp\tn_nonmiss\tmaf\tmac\tqc" > "$CANDS"
echo -e "tag\tcov_mm\tqtl_rank\tdmr_id\tdmr_chr\tdmr_pos\tchosen_cand_rank\tsv_id\tsv_chr\tsv_start\tsv_end\tdist_bp\tused_win_bp\tn_nonmiss\tmaf\tmac\tqc" > "$INCLUDED"
cp "$ID3" "$COVAR_OUT"

: > "$ALL_CAND_IDS"
: > "$ALL_CAND_TSV"

tail -n +2 "$TARGETS" | sort -t$'\t' -k6,6n | \
while IFS=$'\t' read -r meta mm kin idx src_qtl qtl_rank dmr_id chr_raw chr_num pos beta bsd pval; do
  qtl_rank="${qtl_rank//$'\r'/}"
  chr_raw="${chr_raw//$'\r'/}"
  chr_num="${chr_num//$'\r'/}"
  pos="${pos//$'\r'/}"

  if [[ -n "${chr_raw:-}" && "$chr_raw" != "NA" ]]; then
    dmr_chr="$(norm_sl25_chr "$chr_raw")"
  else
    dmr_chr="$(norm_sl25_chr "$chr_num")"
  fi

  cand_rank=0
  while IFS=$'\t' read -r sv_id sv_chr sv_start sv_end dist used_win; do
    [[ -z "${sv_id:-}" ]] && continue
    cand_rank=$((cand_rank+1))
    echo "$sv_id" >> "$ALL_CAND_IDS"
    echo -e "${qtl_rank}\t${dmr_id}\t${dmr_chr}\t${pos}\t${cand_rank}\t${sv_id}\t${sv_chr}\t${sv_start}\t${sv_end}\t${dist}\t${used_win}" >> "$ALL_CAND_TSV"
  done < <(pick_candidates_sv "$SV_COORD" "$dmr_chr" "$pos" "$TOP_VAR_SCAN" || true)
done

sort -u "$ALL_CAND_IDS" > "$TMPDIR/all_candidate_ids.uniq.txt"

if [[ ! -s "$TMPDIR/all_candidate_ids.uniq.txt" ]]; then
  tail -n +2 "$TARGETS" | awk -F'\t' -v OFS='\t' '
    { print "'"$TAG"'", "SV", $6, $7, $8, $10, "NA","NA","NA","NA","NA","NA","NA",0,"NA",0,"FAIL_NO_VARIANT" }' >> "$INCLUDED"
  echo "Wrote:  $COVAR_OUT"
  echo "Chosen: $INCLUDED"
  echo "QC all: $CANDS"
  echo "DONE: $TAG (SV, no candidates)"
  exit 0
fi

bcftools view -i "ID=@${TMPDIR}/all_candidate_ids.uniq.txt" "$SV_VCF" 2>/dev/null \
  | awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $3,$0}' > "$SV_LINES_TSV"

while IFS=$'\t' read -r qtl_rank dmr_id dmr_chr pos cand_rank sv_id sv_chr sv_start sv_end dist used_win; do
  vline="$(get_cached_vline "$sv_id" "$SV_LINES_TSV")"
  [[ -n "$vline" ]] || continue

  dose_col="$TMPDIR/dose.${qtl_rank}.${cand_rank}.txt"
  variant_to_dosage_col "$vline" "$MM_COLMAP" > "$dose_col"

  qc_out="$(qc_dose "$dose_col")"
  qc_flag="$(echo "$qc_out" | awk '{print $1}')"
  n_nonmiss="$(echo "$qc_out" | awk '{print $2}')"
  maf="$(echo "$qc_out" | awk '{print $3}')"
  mac="$(echo "$qc_out" | awk '{print $4}')"

  echo -e "${TAG}\tSV\t${qtl_rank}\t${dmr_id}\t${dmr_chr}\t${pos}\t${cand_rank}\t${sv_id}\t${sv_chr}\t${sv_start}\t${sv_end}\t${dist}\t${used_win}\t${n_nonmiss}\t${maf}\t${mac}\t${qc_flag}" >> "$CANDS"
done < "$ALL_CAND_TSV"

tail -n +2 "$CANDS" | sort -t$'\t' -k3,3n -k12,12n | \
awk -F'\t' -v OFS='\t' '
{
  key=$3 FS $4 FS $5 FS $6
  if(!(key in seen_any)) {
    seen_any[key]=1
    best_any[key]=$0
    best_any_dist[key]=$12+0
    best_any_nonmiss[key]=$14+0
    best_any_mac[key]=$16+0
  } else {
    dist=$12+0; n=$14+0; mac=$16+0
    if(dist < best_any_dist[key] || (dist==best_any_dist[key] && n > best_any_nonmiss[key]) || (dist==best_any_dist[key] && n==best_any_nonmiss[key] && mac > best_any_mac[key])) {
      best_any[key]=$0
      best_any_dist[key]=dist
      best_any_nonmiss[key]=n
      best_any_mac[key]=mac
    }
  }

  if($17=="OK"){
    if(!(key in seen_ok)) {
      seen_ok[key]=1
      best_ok[key]=$0
      best_ok_dist[key]=$12+0
      best_ok_nonmiss[key]=$14+0
      best_ok_mac[key]=$16+0
    } else {
      dist=$12+0; n=$14+0; mac=$16+0
      if(dist < best_ok_dist[key] || (dist==best_ok_dist[key] && n > best_ok_nonmiss[key]) || (dist==best_ok_dist[key] && n==best_ok_nonmiss[key] && mac > best_ok_mac[key])) {
        best_ok[key]=$0
        best_ok_dist[key]=dist
        best_ok_nonmiss[key]=n
        best_ok_mac[key]=mac
      }
    }
  }
}
END{
  for(k in seen_any){
    if(k in seen_ok) print best_ok[k];
    else print best_any[k];
  }
}' | sort -t$'\t' -k3,3n > "$TMPDIR/best_per_dmr.tsv"

while IFS=$'\t' read -r tag cov_mm qtl_rank dmr_id dmr_chr pos cand_rank sv_id sv_chr sv_start sv_end dist used_win n_nonmiss maf mac qc; do
  vline="$(get_cached_vline "$sv_id" "$SV_LINES_TSV")"
  if [[ -n "$vline" ]]; then
    dose_col="$TMPDIR/final_dose.${qtl_rank}.${cand_rank}.txt"
    variant_to_dosage_col "$vline" "$MM_COLMAP" > "$dose_col"
    paste "$COVAR_OUT" "$dose_col" > "$TMPDIR/tmp.cov"
    mv "$TMPDIR/tmp.cov" "$COVAR_OUT"
  elif [[ "$KEEP_FAIL_COLS" == "1" ]]; then
    awk 'BEGIN{OFS="\t"}{print $0,-9}' "$COVAR_OUT" > "$TMPDIR/tmp.cov"
    mv "$TMPDIR/tmp.cov" "$COVAR_OUT"
  fi

  echo -e "${tag}\t${cov_mm}\t${qtl_rank}\t${dmr_id}\t${dmr_chr}\t${pos}\t${cand_rank}\t${sv_id}\t${sv_chr}\t${sv_start}\t${sv_end}\t${dist}\t${used_win}\t${n_nonmiss}\t${maf}\t${mac}\t${qc}" >> "$INCLUDED"
done < "$TMPDIR/best_per_dmr.tsv"

echo "Wrote:  $COVAR_OUT"
echo "Chosen: $INCLUDED"
echo "QC all: $CANDS"
echo "DONE: $TAG (SV)"
