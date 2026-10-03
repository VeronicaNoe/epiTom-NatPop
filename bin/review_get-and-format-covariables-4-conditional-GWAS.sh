#!/usr/bin/env bash
JOB="$(cat "$1")"
TAG="$(echo "$JOB" | cut -d'_' -f1-4)"
COV_MM="$(echo "$JOB" | cut -d'_' -f5)"

WDIR="/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS"
TARGETS="$WDIR/aa_get-target-dmrs/${TAG}.targets.tsv"
TFAM="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/ba_markers/DMR_general_leaf_LD.tfam"

case "$COV_MM" in
  SNP) MM_VCF="/mnt/disk2/vibanez/01_raw-data/01.2_vcfiles/ab_vcf-metabolome/SNP_general_leaf-metabolome_LD.vcf.gz" ;;
  TIP) MM_VCF="$WDIR/tables/TIPs.common.ldpruned.vcf.gz" ;;
  *) echo "ERROR: COV_MM must be SNP or TIP, got: $COV_MM" >&2; exit 2 ;;
esac

TOP_VAR_SCAN="${TOP_VAR_SCAN:-50}"
MIN_NONMISS="${MIN_NONMISS:-50}"
MIN_MAF="${MIN_MAF:-0.01}"
MIN_MAC="${MIN_MAC:-0}"               # set e.g. 5 if you want
KEEP_FAIL_COLS="${KEEP_FAIL_COLS:-0}" # 0=skip failed DMRs (recommended); 1=append -9 for failures

EXPAND_WINS=(50000 100000 250000 500000 1000000 2000000 5000000 10000000)

OUTDIR="$WDIR/ab_covariables/${COV_MM}"
mkdir -p "$OUTDIR"

CANDS="$OUTDIR/${TAG}.cands.tsv"
COVAR_OUT="$OUTDIR/${TAG}__${COV_MM}.covar"
INCLUDED="$OUTDIR/${TAG}__${COV_MM}.included_variants.tsv"

[[ -s "$TARGETS" ]] || { echo "ERROR: missing TARGETS: $TARGETS" >&2; exit 1; }
[[ -s "$TFAM" ]] || { echo "ERROR: missing TFAM: $TFAM" >&2; exit 1; }
[[ -s "$MM_VCF" ]] || { echo "ERROR: missing VCF: $MM_VCF" >&2; exit 1; }
[[ -s "${MM_VCF}.tbi" ]] || { echo "ERROR: missing VCF index: ${MM_VCF}.tbi" >&2; exit 1; }

# ----------------------------
# tfam skeleton + sample order
# ----------------------------
ID3="$OUTDIR/${TAG}.tfam_fid_iid_intercept.tsv"
awk 'BEGIN{OFS="\t"}{print $1,$2,1}' "$TFAM" > "$ID3"

SAMPLES_FILE="$OUTDIR/${TAG}.samples.tfam_order.txt"
awk '{print $2}' "$TFAM" > "$SAMPLES_FILE"

make_colmap() {
  local vcf="$1" out="$2"
  local header
  header="$(zgrep -m1 '^#CHROM' "$vcf")"
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

MM_COLMAP="$OUTDIR/${TAG}.${COV_MM}.colmap.tsv"
make_colmap "$MM_VCF" "$MM_COLMAP"

# ----------------------------
# Chromosome aliases for tabix
# ----------------------------
chrom_aliases() {
  local c="$1"
  # If already non-numeric, also try as-is
  echo "$c"
  # If numeric, try common prefixes
  if [[ "$c" =~ ^[0-9]+$ ]]; then
    printf "chr%s\n" "$c"
    printf "ch%02d\n" "$c"
    printf "SL2.50ch%02d\n" "$c"
    printf "SL2.5ch%02d\n" "$c"
  fi
}

# ----------------------------
# Variant line fetch + dosage + QC
# ----------------------------
get_variant_line_by_refalt() {
  local vcf="$1" chr="$2" pos="$3" ref="$4" alt="$5"
  local lines
  lines="$(tabix "$vcf" "${chr}:${pos}-${pos}" 2>/dev/null || true)"
  [[ -n "$lines" ]] || return 1
  echo "$lines" | awk -v FS="\t" -v ref="$ref" -v alt="$alt" '
    $4==ref && $5==alt {print; found=1; exit}
    END{ if(!found) exit 1 }' && return 0
  echo "$lines" | head -n1
}

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
        split(gtfield,a,":"); g=a[1];
        if(g=="" || g ~ /\./){ print -9; continue; }
        m=split(g,al,/[\/|]/);
        dose=0; bad=0;
        for(j=1;j<=m;j++){
          if(al[j] !~ /^[0-9]+$/){ bad=1; break; }
          if((al[j]+0)!=0) dose+=1;
        }
        if(bad) print -9; else print dose;
      }
    }' /dev/null
}

qc_dose() {
  local dose_file="$1"
  awk -v minN="$MIN_NONMISS" -v minMAF="$MIN_MAF" -v minMAC="$MIN_MAC" '
    $1!=-9 {n++; alt+=$1}
    END{
      if(n==0){print "FAIL 0 NA 0"; exit}
      totA=2*n; ref=totA-alt;
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
# Per-DMR candidate scan (tabix with chr aliases)
# Output: var_chr var_pos var_id ref alt dist used_win
# ----------------------------
pick_candidates() {
  local vcf="$1" chr_num="$2" pos="$3" k="$4"
  local w lo hi chr tmp n

  for w in "${EXPAND_WINS[@]}"; do
    lo=$((pos-w)); ((lo<1)) && lo=1
    hi=$((pos+w))

    for chr in $(chrom_aliases "$chr_num"); do
      tmp="$(mktemp)"
      tabix "$vcf" "${chr}:${lo}-${hi}" 2>/dev/null \
        | awk -v FS="\t" -v OFS="\t" -v pos="$pos" '
            /^#/ {next}
            {
              d = ($2>pos ? $2-pos : pos-$2);
              id=$3; if(id=="." || id=="") id=$1 ":" $2 ":" $4 ":" $5;
              print $1,$2,id,$4,$5,d
            }' \
        | sort -k6,6n | head -n "$k" > "$tmp"
      n=$(wc -l < "$tmp" | tr -d ' ')
      if [[ "$n" -gt 0 ]]; then
        awk -v OFS="\t" -v w="$w" '{print $1,$2,$3,$4,$5,$6,w}' "$tmp"
        rm -f "$tmp"
        return 0
      fi
      rm -f "$tmp"
    done
  done
  return 1
}

# ----------------------------
# Outputs
# ----------------------------
echo -e "tag\tcov_mm\tqtl_rank\tdmr_id\tchr_num\tpos\tcand_rank\tvar_chr\tvar_pos\tvar_id\tref\talt\tdist_bp\tused_win_bp\tn_nonmiss\tmaf\tmac\tqc" > "$CANDS"
echo -e "tag\tcov_mm\tqtl_rank\tdmr_id\tchr_num\tpos\tchosen_cand_rank\tvar_chr\tvar_pos\tvar_id\tref\talt\tdist_bp\tused_win_bp\tn_nonmiss\tmaf\tmac\tqc" > "$INCLUDED"

cp "$ID3" "$COVAR_OUT"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

# Iterate DMRs in qtl_rank order
tail -n +2 "$TARGETS" | sort -t$'\t' -k6,6n | \
while IFS=$'\t' read -r meta mm kin idx src_qtl qtl_rank dmr_id chr_raw chr_num pos beta bsd pval; do
  qtl_rank="${qtl_rank//$'\r'/}"
  chr_num="${chr_num//$'\r'/}"
  pos="${pos//$'\r'/}"

  mapfile -t cand_lines < <(pick_candidates "$MM_VCF" "$chr_num" "$pos" "$TOP_VAR_SCAN" || true)

  best_ok_set=0
  best_any_set=0

  # best OK
  best_ok_nonmiss=-1; best_ok_mac=-1; best_ok_dist=999999999
  best_ok_cand_rank="NA"; best_ok_line=""; best_ok_dose=""; best_ok_maf="NA"

  # best ANY (even FAIL)
  best_any_nonmiss=-1; best_any_mac=-1; best_any_dist=999999999
  best_any_cand_rank="NA"; best_any_line=""; best_any_dose=""; best_any_qc="FAIL"; best_any_maf="NA"

  cand_rank=0
  for cl in "${cand_lines[@]}"; do
    cand_rank=$((cand_rank+1))
    IFS=$'\t' read -r var_chr var_pos var_id ref alt dist used_win <<< "$cl"

    vline="$(get_variant_line_by_refalt "$MM_VCF" "$var_chr" "$var_pos" "$ref" "$alt" || true)"
    [[ -n "$vline" ]] || continue

    dose_col="$TMPDIR/dose.${qtl_rank}.${cand_rank}.txt"
    variant_to_dosage_col "$vline" "$MM_COLMAP" > "$dose_col"

    qc_out="$(qc_dose "$dose_col")"
    qc_flag="$(echo "$qc_out" | awk '{print $1}')"
    n_nonmiss="$(echo "$qc_out" | awk '{print $2}')"
    maf="$(echo "$qc_out" | awk '{print $3}')"
    mac="$(echo "$qc_out" | awk '{print $4}')"

    echo -e "${TAG}\t${COV_MM}\t${qtl_rank}\t${dmr_id}\t${chr_num}\t${pos}\t${cand_rank}\t${var_chr}\t${var_pos}\t${var_id}\t${ref}\t${alt}\t${dist}\t${used_win}\t${n_nonmiss}\t${maf}\t${mac}\t${qc_flag}" >> "$CANDS"

    # best ANY: max nonmiss, then max mac, then min dist
    if (( n_nonmiss > best_any_nonmiss )) || \
       (( n_nonmiss == best_any_nonmiss && mac > best_any_mac )) || \
       (( n_nonmiss == best_any_nonmiss && mac == best_any_mac && dist < best_any_dist )); then
      best_any_set=1
      best_any_nonmiss=$n_nonmiss
      best_any_mac=$mac
      best_any_dist=$dist
      best_any_maf=$maf
      best_any_qc=$qc_flag
      best_any_cand_rank=$cand_rank
      best_any_line="$cl"
      best_any_dose="$dose_col"
    fi

    # best OK only among QC-pass
    if [[ "$qc_flag" == "OK" ]]; then
      if (( n_nonmiss > best_ok_nonmiss )) || \
         (( n_nonmiss == best_ok_nonmiss && mac > best_ok_mac )) || \
         (( n_nonmiss == best_ok_nonmiss && mac == best_ok_mac && dist < best_ok_dist )); then
        best_ok_set=1
        best_ok_nonmiss=$n_nonmiss
        best_ok_mac=$mac
        best_ok_dist=$dist
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

    IFS=$'\t' read -r var_chr var_pos var_id ref alt dist used_win <<< "$chosen_line"
    echo -e "${TAG}\t${COV_MM}\t${qtl_rank}\t${dmr_id}\t${chr_num}\t${pos}\t${chosen_rank}\t${var_chr}\t${var_pos}\t${var_id}\t${ref}\t${alt}\t${dist}\t${used_win}\t${chosen_nonmiss}\t${chosen_maf}\t${chosen_mac}\t${chosen_qc}" >> "$INCLUDED"
  else
    # no usable variant → either skip column or append -9 column
    if [[ "$KEEP_FAIL_COLS" == "1" ]]; then
      awk 'BEGIN{OFS="\t"}{print $0,-9}' "$COVAR_OUT" > "$TMPDIR/tmp.cov"
      mv "$TMPDIR/tmp.cov" "$COVAR_OUT"
    fi
    echo -e "${TAG}\t${COV_MM}\t${qtl_rank}\t${dmr_id}\t${chr_num}\t${pos}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0\tNA\t0\t${chosen_qc}" >> "$INCLUDED"
  fi
done

echo "Wrote:  $COVAR_OUT"
echo "Chosen: $INCLUDED"
echo "QC all: $CANDS"
echo "DONE: $TAG ($COV_MM)"
