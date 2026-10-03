#!/usr/bin/env bash
set -euo pipefail

############################################
# CONFIG
############################################
SIG_GWAS="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results/sig"
WDIR="/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS"
COVAR_DIR="$WDIR/aa_covariables"
BASE_COVAR_MODE="per_meta"
GLOBAL_COVAR=""

OUTDIR="$WDIR/ac_outputs"
mkdir -p "$OUTDIR" "$COVAR_DIR" "$WDIR/tables" "$WDIR/logs"

# Variant sources (must be bgzip+tabix indexed .vcf.gz)
SNP_VCF="/mnt/disk2/vibanez/01_raw-data/01.2_vcfiles/ab_vcf-metabolome/SNP_general_leaf-metabolome_LD.vcf.gz"
TIP_VCF="$WDIR/tables/TIPs.common.ldpruned.vcf.gz"

# Lead selection + conditioning parameters
PTHRESH="5e-8"     # used to pick "lead DMRs"; if none pass, falls back to global top N
TOP_DMR_N=1        # how many lead DMRs per metabolite
MAX_WIN_BP=10000000  # +/- window for nearest variants (10 Mb)
TOP_VAR_K=3        # how many nearest variants per source (SNP and TIP)
EXPAND_WINS=(50000 100000 250000 500000 1000000 2000000 5000000 10000000)  # progressive search
############################################
# HELPERS
############################################
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command: $1" >&2; exit 1; }; }
check_vcf_index() {
  local v="$1"
  [[ -f "$v" ]] || { echo "ERROR: VCF not found: $v" >&2; exit 1; }
  [[ "$v" =~ \.vcf\.gz$ ]] || { echo "ERROR: VCF must be bgzipped .vcf.gz for tabix: $v" >&2; exit 1; }
  [[ -f "${v}.tbi" ]] || { echo "ERROR: missing tabix index ${v}.tbi (run: tabix -p vcf $v)" >&2; exit 1; }
}
# 1) bgzip-compress (creates .vcf.gz)
#bgzip -c SNP_general_leaf-metabolome_LD.vcf > SNP_general_leaf-metabolome_LD.vcf.gz
# 2) index for random access
#tabix -p vcf SNP_general_leaf-metabolome_LD.vcf.gz
# quick check
#tabix -H SNP_general_leaf-metabolome_LD.vcf.gz | head
#bgzip -c TIPs.vcf > TIPs.vcf.gz
# 2) index for random access
#tabix -p vcf TIPs.vcf.gz
# quick check
#tabix -H TIPs.vcf.gz | head
# Convert DMR chrom string to VCF chrom (your SNP VCF uses numeric chrom like "1")
# Examples:
#   ch01 -> 1
#   ch00 -> 0
#   chr1 -> 1
chr_to_num() {
  local c="$1"
  c="${c#chr}"
  c="${c#ch}"
  # strip leading zeros but keep "0"
  c="$(echo "$c" | sed 's/^0\+\([0-9]\)/\1/; s/^0$/0/')"
  echo "$c"
}

# Make a sample order list from the DMR tfam (IID column)
SAMPLES_FILE="$WDIR/tables/samples.tfam_order.txt"
awk '{print $2}' "/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/ba_markers/DMR_general_leaf_LD.tfam" > "$SAMPLES_FILE"

# Build a "sample -> VCF column index" map restricted to our sample set (for fast dosage extraction)
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
    { s=$1; if(s in col) print s, col[s]; else print s, 0; }
  ' "$SAMPLES_FILE" > "$out"
}

# Extract variant line at chr:pos, and optionally match ID (if multiple at same pos)
get_variant_line() {
  local vcf="$1" chr="$2" pos="$3" id="$4"
  if [[ -z "$id" || "$id" == "." ]]; then
    tabix "$vcf" "${chr}:${pos}-${pos}" | head -n1
  else
    tabix "$vcf" "${chr}:${pos}-${pos}" | awk -v id="$id" 'BEGIN{FS="\t"} $3==id {print; exit}'
  fi
}

# Convert a single VCF variant line to a dosage column in tfam order using colmap
# Dosage = count of non-zero alleles in GT (0/0=0, 0/1=1, 1/1=2; missing -> NA)
variant_to_dosage_col() {
  local variant_line="$1"
  local colmap="$2"
  awk -v line="$variant_line" -v colmap="$colmap" '
    BEGIN{
      FS=OFS="\t";
      while((getline<colmap)>0){
        samp[++n]=$1; idx[$1]=$2;
      }
      close(colmap);
      split(line, F, "\t");
      for(i=1;i<=n;i++){
        c=idx[samp[i]];
        if(c==0){ print "NA"; continue; }
        gt=F[c];
        split(gt, a, ":"); g=a[1];
        if(g ~ /\./){ print "NA"; continue; }
        m=split(g, al, /[\/|]/);
        dose=0;
        bad=0;
        for(j=1;j<=m;j++){
          if(al[j]=="."){ bad=1; break; }
          if(al[j]!="0") dose+=1;
        }
        if(bad) print "NA"; else print dose;
      }
    }
  ' /dev/null
}

# Pick K nearest variants to a position, searching progressively up to MAX_WIN_BP
# Output: type  meta  lead_dmr  vcf_chr  vcf_pos  vcf_id  dist_bp  used_win_bp
pick_nearest() {
  local vcf="$1" vtype="$2" meta="$3" lead="$4" chr="$5" pos="$6" k="$7"
  local w lo hi out tmp
  tmp="$(mktemp)"
  for w in "${EXPAND_WINS[@]}"; do
    lo=$((pos-w)); ((lo<1)) && lo=1
    hi=$((pos+w))
    # stream region; keep top-k by distance without sorting entire list
    tabix "$vcf" "${chr}:${lo}-${hi}" 2>/dev/null \
      | awk -v FS="\t" -v OFS="\t" -v pos="$pos" -v k="$k" -v w="$w" -v meta="$meta" -v lead="$lead" -v vtype="$vtype" '
          BEGIN{n=0}
          /^#/ {next}
          {
            d = ($2>pos ? $2-pos : pos-$2);
            id=$3; if(id=="." || id=="") id=$1 ":" $2 ":" $4 ":" $5;

            # insert into best list (size <= k)
            if(n < k){
              n++; bestd[n]=d; bestchr[n]=$1; bestpos[n]=$2; bestid[n]=id;
            } else {
              # find current worst
              wi=1; wd=bestd[1];
              for(i=2;i<=n;i++) if(bestd[i]>wd){wd=bestd[i]; wi=i;}
              if(d < wd){
                bestd[wi]=d; bestchr[wi]=$1; bestpos[wi]=$2; bestid[wi]=id;
              }
            }
          }
          END{
            for(i=1;i<=n;i++) print vtype, meta, lead, bestchr[i], bestpos[i], bestid[i], bestd[i], w;
          }
        ' >> "$tmp"

    # if we got at least k lines, stop expanding
    if [[ "$(wc -l < "$tmp")" -ge "$k" ]]; then break; fi
  done
  cat "$tmp"
  rm -f "$tmp"
}

# Extract baseline stats for a specific DMR from the baseline .ps.gz
baseline_stats() {
  local ps_gz="$1" dmr="$2"
  zcat "$ps_gz" | awk -v d="$dmr" 'BEGIN{FS="\t"} $1==d {print $2"\t"$3"\t"$4; exit}'
}

# Extract conditional stats for a specific DMR from an emmax .ps file (uncompressed)
cond_stats() {
  local ps="$1" dmr="$2"
  awk -v d="$dmr" 'BEGIN{FS="\t"} $1==d {print $2"\t"$3"\t"$4; exit}' "$ps"
}

############################################
# PRECHECKS
############################################
need_cmd awk
need_cmd zgrep
need_cmd zcat
need_cmd tabix
need_cmd pigz

check_vcf_index "$SNP_VCF"
check_vcf_index "$TIP_VCF"

SNP_COLMAP="$WDIR/tables/snp.colmap.tsv"
TIP_COLMAP="$WDIR/tables/tip.colmap.tsv"
make_colmap "$SNP_VCF" "$SNP_COLMAP"
make_colmap "$TIP_VCF" "$TIP_COLMAP"
