#!/usr/bin/env bash
TAG="$(cat "$1")"
META="$(cat "$1" | cut -d'.' -f1)"
INDEX="$(cat "$1" | cut -d'_' -f4)"
MM="$(cat "$1" | cut -d'_' -f2)"
KINSHIP="$(cat "$1" | cut -d'_' -f3)"

SIG_DIR="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results/sig"
WDIR="/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS"
TDIR="$WDIR/aa_get-target-dmrs"

TARGETS="$TDIR/ba_targets/${TAG}.target"
TFAM="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/ba_markers/DMR_general_leaf_LD.tfam"

SV_VCF="/mnt/disk2/vibanez/otherAnalysis/13_DMR-GWAS-SV/aa_markers/graphpanSV.sharedSamples.vcf.gz"
SV_COORD="/mnt/disk2/vibanez/otherAnalysis/review_DMR-GWAS/aa_DMR-GWAS-SV_SL25/sv_id2coord_sl25_primary.tsv"

TOP_VAR_SCAN="${TOP_VAR_SCAN:-50}"
MIN_NONMISS="${MIN_NONMISS:-50}"
MIN_MAF="${MIN_MAF:-0.01}"
MIN_MAC="${MIN_MAC:-0}"
KEEP_FAIL_COLS="${KEEP_FAIL_COLS:-0}"

EXPAND_WINS=(50000 100000 250000 500000 1000000 2000000 5000000 10000000)

OUTDIR="$WDIR/ab_covariables/SV"
mkdir -p "$OUTDIR"

CANDS="$OUTDIR/${TAG}.cands.tsv"
COVAR_OUT="$OUTDIR/${TAG}__SV.covar"
INCLUDED="$OUTDIR/${TAG}__SV.included_variants.tsv"

[[ -s "$TARGETS" ]] || { echo "ERROR: missing TARGETS: $TARGETS" >&2; exit 1; }
[[ -s "$TFAM" ]] || { echo "ERROR: missing TFAM: $TFAM" >&2; exit 1; }
[[ -s "$SV_VCF" ]] || { echo "ERROR: missing SV_VCF: $SV_VCF" >&2; exit 1; }
[[ -s "${SV_VCF}.tbi" ]] || { echo "ERROR: missing VCF index: ${SV_VCF}.tbi" >&2; exit 1; }
[[ -s "$SV_COORD" ]] || { echo "ERROR: missing SV_COORD: $SV_COORD" >&2; exit 1; }

# ----------------------------
# tfam skeleton + sample order
# ----------------------------
ID3="$OUTDIR/${TAG}.tfam_fid_iid_intercept.tsv"
awk 'BEGIN{OFS="\t"}{print $1,$2,1}' "$TFAM" > "$ID3"

SAMPLES_FILE="$OUTDIR/${TAG}.samples.tfam_order.txt"
awk '{print $2}' "$TFAM" > "$SAMPLES_FILE"

# ----------------------------
# Build sample -> VCF column map
# If sample not in VCF, assign 0 and dosage becomes -9
# ----------------------------
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

MM_COLMAP="$OUTDIR/${TAG}.SV.colmap.tsv"
make_colmap "$SV_VCF" "$MM_COLMAP"

# ----------------------------
# Chromosome normalization to SL2.50chXX
# ----------------------------
norm_sl25_chr() {
  local raw="$1"
  raw="${raw//$'\r'/}"

  case "$raw" in
    SL2.50ch*) echo "$raw" ;;
    SL2.5ch*)
      echo "${raw/SL2.5ch/SL2.50ch}"
      ;;
    ch[0-9][0-9])
      printf "SL2.50%s\n" "$raw"
      ;;
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

# ----------------------------
# Fetch VCF line by ID
# ----------------------------
get_variant_line_by_id() {
  local vcf="$1" id="$2"
  bcftools view -i "ID=\"${id}\"" "$vcf" 2>/dev/null | awk 'BEGIN{FS="\t"} !/^#/ {print; exit}'
}

# ----------------------------
# Convert GT to dosage in TFAM order
# non-zero allele count => dosage
# missing or absent sample => -9
# ----------------------------
variant_to_dosage_col() {
  local variant_line="$1" colmap="$2"
  awk -v line="$variant_line" -v colmap="$colmap" '
    BEGIN{
      FS=OFS="\t";
      while((getline<colmap)>0){ samp[++n]=$1; idx[$1]=$2; }
      close(colmap);
      split(line,F,"\t");

      for(i=1;i<=n;i++){
        c=idx[samp[i]];
        if(c==0){ print -9; continue; }

        gtfield=F[c];
        split(gtfield,a,":");
        g=a[1];

        if(g=="" || g ~ /\./){ print -9; continue; }

        m=split(g,al,/[\/|]/);
        dose=0; bad=0;
        for(j=1;j<=m;j++){
          if(al[j] !~ /^[0-9]+$/){ bad=1; break; }
          if((al[j]+0)!=0) dose+=1;
        }

        if(bad) print -9;
        else    print dose;
      }
    }' /dev/null
}

# ----------------------------
# QC dosage column
# ----------------------------
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

# ----------------------------
# Pick candidate SVs in SL2.5 space
# Input coordinate file columns assumed:
# 1=sv_id  2=chr_num  3=center  4=chr_sl25  5=start_sl25  6=end_sl25
#
# Output:
# sv_id  sv_chr  sv_start  sv_end  dist  used_win
# ----------------------------
pick_candidates_sv() {
  local coord_tsv="$1" dmr_chr="$2" dmr_pos="$3" k="$4"
  local w tmp n

  for w in "${EXPAND_WINS[@]}"; do
    tmp="$(mktemp)"
    awk -v FS="\t" -v OFS="\t" -v chr="$dmr_chr" -v pos="$dmr_pos" -v w="$w" '
      $4==chr {
        sv_id=$1;
        sv_chr=$4;
        sv_start=$5+0;
        sv_end=$6+0;

        if(pos < sv_start)      dist = sv_start - pos;
        else if(pos > sv_end)   dist = pos - sv_end;
        else                    dist = 0;

        if(dist <= w){
          print sv_id, sv_chr, sv_start, sv_end, dist;
        }
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

# ----------------------------
# Outputs
# ----------------------------
echo -e "tag\tcov_mm\tqtl_rank\tdmr_id\tdmr_chr\tdmr_pos\tcand_rank\tsv_id\tsv_chr\tsv_start\tsv_end\tdist_bp\tused_win_bp\tn_nonmiss\tmaf\tmac\tqc" > "$CANDS"
echo -e "tag\tcov_mm\tqtl_rank\tdmr_id\tdmr_chr\tdmr_pos\tchosen_cand_rank\tsv_id\tsv_chr\tsv_start\tsv_end\tdist_bp\tused_win_bp\tn_nonmiss\tmaf\tmac\tqc" > "$INCLUDED"

cp "$ID3" "$COVAR_OUT"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

# ----------------------------
# Iterate DMRs in qtl_rank order
# Expected TARGETS columns as in your old script:
# meta mm kin idx src_qtl qtl_rank dmr_id chr_raw chr_num pos beta bsd pval
# ----------------------------
tail -n +2 "$TARGETS" | sort -t$'\t' -k6,6n | \
while IFS=$'\t' read -r meta mm kin idx src_qtl qtl_rank dmr_id chr_raw chr_num pos beta bsd pval; do
  qtl_rank="${qtl_rank//$'\r'/}"
  chr_raw="${chr_raw//$'\r'/}"
  chr_num="${chr_num//$'\r'/}"
  pos="${pos//$'\r'/}"

  # Prefer chr_raw if already informative, otherwise chr_num
  if [[ -n "${chr_raw:-}" && "$chr_raw" != "NA" ]]; then
    dmr_chr="$(norm_sl25_chr "$chr_raw")"
  else
    dmr_chr="$(norm_sl25_chr "$chr_num")"
  fi

  mapfile -t cand_lines < <(pick_candidates_sv "$SV_COORD" "$dmr_chr" "$pos" "$TOP_VAR_SCAN" || true)

  best_ok_set=0
  best_any_set=0

  # best OK: prioritize distance, then n_nonmiss, then MAC
  best_ok_dist=999999999
  best_ok_nonmiss=-1
  best_ok_mac=-1
  best_ok_cand_rank="NA"
  best_ok_line=""
  best_ok_dose=""
  best_ok_maf="NA"

  # best ANY (even FAIL): same ranking
  best_any_dist=999999999
  best_any_nonmiss=-1
  best_any_mac=-1
  best_any_cand_rank="NA"
  best_any_line=""
  best_any_dose=""
  best_any_qc="FAIL"
  best_any_maf="NA"

  cand_rank=0
  for cl in "${cand_lines[@]}"; do
    cand_rank=$((cand_rank+1))
    IFS=$'\t' read -r sv_id sv_chr sv_start sv_end dist used_win <<< "$cl"

    vline="$(get_variant_line_by_id "$SV_VCF" "$sv_id" || true)"
    [[ -n "$vline" ]] || continue

    dose_col="$TMPDIR/dose.${qtl_rank}.${cand_rank}.txt"
    variant_to_dosage_col "$vline" "$MM_COLMAP" > "$dose_col"

    qc_out="$(qc_dose "$dose_col")"
    qc_flag="$(echo "$qc_out" | awk '{print $1}')"
    n_nonmiss="$(echo "$qc_out" | awk '{print $2}')"
    maf="$(echo "$qc_out" | awk '{print $3}')"
    mac="$(echo "$qc_out" | awk '{print $4}')"

    echo -e "${TAG}\tSV\t${qtl_rank}\t${dmr_id}\t${dmr_chr}\t${pos}\t${cand_rank}\t${sv_id}\t${sv_chr}\t${sv_start}\t${sv_end}\t${dist}\t${used_win}\t${n_nonmiss}\t${maf}\t${mac}\t${qc_flag}" >> "$CANDS"

    # best ANY
    if (( dist < best_any_dist )) || \
       (( dist == best_any_dist && n_nonmiss > best_any_nonmiss )) || \
       (( dist == best_any_dist && n_nonmiss == best_any_nonmiss && mac > best_any_mac )); then
      best_any_set=1
      best_any_dist=$dist
      best_any_nonmiss=$n_nonmiss
      best_any_mac=$mac
      best_any_maf=$maf
      best_any_qc=$qc_flag
      best_any_cand_rank=$cand_rank
      best_any_line="$cl"
      best_any_dose="$dose_col"
    fi

    # best OK
    if [[ "$qc_flag" == "OK" ]]; then
      if (( dist < best_ok_dist )) || \
         (( dist == best_ok_dist && n_nonmiss > best_ok_nonmiss )) || \
         (( dist == best_ok_dist && n_nonmiss == best_ok_nonmiss && mac > best_ok_mac )); then
        best_ok_set=1
        best_ok_dist=$dist
        best_ok_nonmiss=$n_nonmiss
        best_ok_mac=$mac
        best_ok_maf=$maf
        best_ok_cand_rank=$cand_rank
        best_ok_line="$cl"
        best_ok_dose="$dose_col"
      fi
    fi
  done

  # choose best OK if exists; else best ANY if exists
  if (( best_ok_set == 1 )); then
    chosen_rank="$best_ok_cand_rank"
    chosen_line="$best_ok_line"
    chosen_dose="$best_ok_dose"
    chosen_nonmiss="$best_ok_nonmiss"
    chosen_maf="$best_ok_maf"
    chosen_mac="$best_ok_mac"
    chosen_qc="OK"
  elif (( best_any_set == 1 )); then
    chosen_rank="$best_any_cand_rank"
    chosen_line="$best_any_line"
    chosen_dose="$best_any_dose"
    chosen_nonmiss="$best_any_nonmiss"
    chosen_maf="$best_any_maf"
    chosen_mac="$best_any_mac"
    chosen_qc="$best_any_qc"
  else
    chosen_rank="NA"
    chosen_line=""
    chosen_dose=""
    chosen_nonmiss=0
    chosen_maf="NA"
    chosen_mac=0
    chosen_qc="FAIL_NO_VARIANT"
  fi

  if [[ -n "$chosen_dose" ]]; then
    paste "$COVAR_OUT" "$chosen_dose" > "$TMPDIR/tmp.cov"
    mv "$TMPDIR/tmp.cov" "$COVAR_OUT"

    IFS=$'\t' read -r sv_id sv_chr sv_start sv_end dist used_win <<< "$chosen_line"
    echo -e "${TAG}\tSV\t${qtl_rank}\t${dmr_id}\t${dmr_chr}\t${pos}\t${chosen_rank}\t${sv_id}\t${sv_chr}\t${sv_start}\t${sv_end}\t${dist}\t${used_win}\t${chosen_nonmiss}\t${chosen_maf}\t${chosen_mac}\t${chosen_qc}" >> "$INCLUDED"
  else
    if [[ "$KEEP_FAIL_COLS" == "1" ]]; then
      awk 'BEGIN{OFS="\t"}{print $0,-9}' "$COVAR_OUT" > "$TMPDIR/tmp.cov"
      mv "$TMPDIR/tmp.cov" "$COVAR_OUT"
    fi
    echo -e "${TAG}\tSV\t${qtl_rank}\t${dmr_id}\t${dmr_chr}\t${pos}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0\tNA\t0\t${chosen_qc}" >> "$INCLUDED"
  fi
done
