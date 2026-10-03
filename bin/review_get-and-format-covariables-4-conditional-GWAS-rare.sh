#!/usr/bin/env bash
JOB="$(cat "$1")"
TAG="$(echo "$JOB" | cut -d'_' -f1-4)"
COV_MM="$(echo "$JOB" | cut -d'_' -f5)"
WDIR="/mnt/disk2/vibanez/otherAnalysis/09_conditional-metabolite-GWAS-rare-variants"
TARGETS="$WDIR/aa_get-target-dmrs/${TAG}.targets.tsv"
TFAM="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/ba_markers/DMR_general_leaf_LD.tfam"

case "$COV_MM" in
  SNP) MM_VCF="$WDIR/tables/SNP_rare_LDpruned.vcf.gz" ;;
  TIP) MM_VCF="/mnt/disk2/vibanez/otherAnalysis/02_TIPs/TIPs.vcf.gz" ;;
  *) echo "ERROR: COV_MM must be SNP or TIP, got: $COV_MM" >&2; exit 2 ;;
esac

TOP_VAR_SCAN=50
TOP_VAR_KEEP=1

# Keep rare variants: allow very small MAF, but avoid pathological singletons with MIN_MAC
MIN_NONMISS=50
MIN_MAF="${MIN_MAF:-0.0}"     # set to 0.001 if you want minimal filtering
MIN_MAC="${MIN_MAC:-3}"       # minor allele count among non-missing (e.g. 3 or 5)

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

# tfam skeleton: FID IID 1
ID3="$OUTDIR/${TAG}.tfam_fid_iid_intercept.tsv"
awk 'BEGIN{OFS="\t"}{print $1,$2,1}' "$TFAM" > "$ID3"

SAMPLES_FILE="$OUTDIR/${TAG}.samples.tfam_order.txt"
awk '{print $2}' "$TFAM" > "$SAMPLES_FILE"

make_colmap() {
  local vcf="$1"
  local out="$2"
  local header
  header="$(zgrep -m1 '^#CHROM' "$vcf")"
  awk -v hdr="$header" '
    BEGIN{
      FS=OFS="\t";
      n=split(hdr, H, "\t");
      for(i=1;i<=n;i++) col[H[i]]=i;
    }
    function strip_suffix(s){ sub(/_[0-9]+$/, "", s); return s; }
    {
      s=$1;
      if(s in col) { print s, col[s]; next; }
      s2=strip_suffix(s);
      if(s2 in col) { print $1, col[s2]; next; }
      print $1, 0;
    }
  ' "$SAMPLES_FILE" > "$out"
}

MM_COLMAP="$OUTDIR/${TAG}.${COV_MM}.colmap.tsv"
make_colmap "$MM_VCF" "$MM_COLMAP"

pick_nearest() {
  local vcf="$1" vtype="$2" meta="$3" lead="$4" chr="$5" pos="$6" k="$7"
  local w lo hi tmp
  tmp="$(mktemp)"
  for w in "${EXPAND_WINS[@]}"; do
    lo=$((pos-w)); ((lo<1)) && lo=1
    hi=$((pos+w))
    tabix "$vcf" "${chr}:${lo}-${hi}" 2>/dev/null \
      | awk -v FS="\t" -v OFS="\t" -v pos="$pos" -v k="$k" -v w="$w" -v meta="$meta" -v lead="$lead" -v vtype="$vtype" '
          BEGIN{n=0}
          /^#/ {next}
          {
            d = ($2>pos ? $2-pos : pos-$2);
            id=$3; if(id=="." || id=="") id=$1 ":" $2 ":" $4 ":" $5;
            ref=$4; alt=$5;

            if(n < k){
              n++; bestd[n]=d; bestchr[n]=$1; bestpos[n]=$2; bestid[n]=id; bestref[n]=ref; bestalt[n]=alt;
            } else {
              wi=1; wd=bestd[1];
              for(i=2;i<=n;i++) if(bestd[i]>wd){wd=bestd[i]; wi=i;}
              if(d < wd){
                bestd[wi]=d; bestchr[wi]=$1; bestpos[wi]=$2; bestid[wi]=id; bestref[wi]=ref; bestalt[wi]=alt;
              }
            }
          }
          END{
            for(i=1;i<=n;i++) print vtype, meta, lead, bestchr[i], bestpos[i], bestid[i], bestref[i], bestalt[i], bestd[i], w;
          }
        ' >> "$tmp"
    [[ "$(wc -l < "$tmp")" -ge "$k" ]] && break
  done
  cat "$tmp"
  rm -f "$tmp"
}

get_variant_line_by_refalt() {
  local vcf="$1" chr="$2" pos="$3" ref="$4" alt="$5"
  local lines
  lines="$(tabix "$vcf" "${chr}:${pos}-${pos}" 2>/dev/null || true)"
  [[ -n "$lines" ]] || return 1
  echo "$lines" | awk -v FS="\t" -v ref="$ref" -v alt="$alt" '
    $4==ref && $5==alt {print; found=1; exit}
    END{ if(!found) exit 1 }
  ' && return 0
  echo "$lines" | head -n1
}

variant_to_dosage_col() {
  local variant_line="$1"
  local colmap="$2"
  awk -v line="$variant_line" -v colmap="$colmap" '
    BEGIN{
      FS=OFS="\t";
      while((getline<colmap)>0){ samp[++n]=$1; idx[$1]=$2; }
      close(colmap);
      split(line, F, "\t");
      for(i=1;i<=n;i++){
        c=idx[samp[i]];
        if(c==0){ print -9; continue; }
        gtfield=F[c];
        split(gtfield, a, ":"); g=a[1];
        if(g=="" || g ~ /\./){ print -9; continue; }
        m=split(g, al, /[\/|]/);
        dose=0; bad=0;
        for(j=1;j<=m;j++){
          if(al[j] !~ /^[0-9]+$/){ bad=1; break; }
          if((al[j]+0)!=0) dose+=1;
        }
        if(bad) print -9; else print dose;
      }
    }
  ' /dev/null
}

# read lead(s) from targets: meta, lead_dmr, chr_num, pos
mapfile -t LEADS < <(tail -n +2 "$TARGETS" | awk -F'\t' 'NF>=9{print $1"\t"$6"\t"$8"\t"$9}')

echo -e "tag\tcov_mm\tmeta\tlead_dmr\tvar_chr\tvar_pos\tvar_id\tref\talt\tdist_bp\tused_win_bp\trank" > "$CANDS"

for line in "${LEADS[@]}"; do
  meta="$(echo "$line" | cut -f1)"
  lead="$(echo "$line" | cut -f2)"
  chr_num="$(echo "$line" | cut -f3)"
  pos="$(echo "$line" | cut -f4)"

  tmp_all="$(mktemp)"
  pick_nearest "$MM_VCF" "$COV_MM" "$meta" "$lead" "$chr_num" "$pos" "$TOP_VAR_SCAN" >> "$tmp_all" || true
  tmp_sorted="$(mktemp)"
  awk 'BEGIN{FS="\t"} NF>=10{print}' "$tmp_all" | sort -k9,9g > "$tmp_sorted"

  awk -v tag="$TAG" -v cov="$COV_MM" -v OFS="\t" '
    BEGIN{FS="\t"; r=0}
    { r++; print tag,cov,$2,$3,$4,$5,$6,$7,$8,$9,$10,r; }
  ' "$tmp_sorted" >> "$CANDS"

  rm -f "$tmp_all" "$tmp_sorted"
done

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

cp "$ID3" "$COVAR_OUT"
echo -e "tag\tcov_mm\tkept_rank\tcand_rank\tvar_chr\tvar_pos\tvar_id\tref\talt\tdist_bp\tused_win_bp\tn_nonmiss\tmaf\tmac" > "$INCLUDED"

kept=0

while IFS=$'\t' read -r tag cov_mm meta lead var_chr var_pos var_id ref alt dist used_win cand_rank; do
  [[ $kept -ge $TOP_VAR_KEEP ]] && break

  vline="$(get_variant_line_by_refalt "$MM_VCF" "$var_chr" "$var_pos" "$ref" "$alt" || true)"
  [[ -n "$vline" ]] || continue

  dose_col="$TMPDIR/dose.cand${cand_rank}.txt"
  variant_to_dosage_col "$vline" "$MM_COLMAP" > "$dose_col"

  qc="$(awk -v minN="$MIN_NONMISS" -v minMAF="$MIN_MAF" -v minMAC="$MIN_MAC" '
    $1!=-9 {n++; s+=$1}
    END{
      if(n==0){print "FAIL\t0\tNA\t0"; exit}
      p=(s/(2*n));
      maf=(p<0.5? p : 1-p);
      mac=(maf*2*n);
      if(n<minN || maf<minMAF || mac<minMAC){printf("FAIL\t%d\t%.6g\t%.0f\n", n, maf, mac); exit}
      printf("OK\t%d\t%.6g\t%.0f\n", n, maf, mac)
    }' "$dose_col")"

  status="$(echo "$qc" | cut -f1)"
  n_nonmiss="$(echo "$qc" | cut -f2)"
  maf="$(echo "$qc" | cut -f3)"
  mac="$(echo "$qc" | cut -f4)"

  [[ "$status" == "OK" ]] || continue

  paste "$COVAR_OUT" "$dose_col" > "$TMPDIR/tmp.cov"
  mv "$TMPDIR/tmp.cov" "$COVAR_OUT"

  kept=$((kept+1))
  echo -e "${TAG}\t${COV_MM}\t${kept}\t${cand_rank}\t${var_chr}\t${var_pos}\t${var_id}\t${ref}\t${alt}\t${dist}\t${used_win}\t${n_nonmiss}\t${maf}\t${mac}" >> "$INCLUDED"
done < <(tail -n +2 "$CANDS" | sort -t$'\t' -k12,12n)

# pad with -9 columns if nothing passed QC
while [[ $kept -lt $TOP_VAR_KEEP ]]; do
  awk 'BEGIN{OFS="\t"}{print $0,-9}' "$COVAR_OUT" > "$TMPDIR/tmp.cov"
  mv "$TMPDIR/tmp.cov" "$COVAR_OUT"
  kept=$((kept+1))
done

echo "Wrote:      $COVAR_OUT"
echo "Included:   $INCLUDED"
echo "Candidates: $CANDS"
echo "DONE: $TAG ($COV_MM)"
