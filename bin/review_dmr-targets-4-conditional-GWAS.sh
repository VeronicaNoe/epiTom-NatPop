#!/usr/bin/env bash
set -euo pipefail

TAG="$(cat "$1")"
META="$(cat "$1" | cut -d'.' -f1)"
INDEX="$(cat "$1" | cut -d'_' -f4)"
MM="$(cat "$1" | cut -d'_' -f2)"
KINSHIP="$(cat "$1" | cut -d'_' -f3)"

SIG_DIR="/mnt/disk2/vibanez/10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results/sig"
WDIR="/mnt/disk2/vibanez/otherAnalysis/06_conditional-metabolite-GWAS"
TDIR="$WDIR/aa_get-target-dmrs"
mkdir -p "$TDIR"

# cap (0 = keep all)
MAX_QTL="${MAX_QTL:-0}"

QTL_FILE="$SIG_DIR/${TAG}.QTL"
OUT="$TDIR/${TAG}.targets.tsv"

[[ -s "$QTL_FILE" ]] || { echo "ERROR: missing QTL file: $QTL_FILE" >&2; exit 1; }

echo -e "meta\tmm\tkinship_label\tindex\tsrc_qtl\tqtl_rank\tdmr_id\tchr_raw\tchr_num\tpos\tbeta\tbeta_sd\tp" > "$OUT"

cat "$QTL_FILE" \
  | awk 'BEGIN{FS="[ \t]+"; OFS="\t"} NF>=4 {print $1,$2,$3,$4}' \
  | sort -gk4,4 \
  | awk -v OFS="\t" -v meta="$META" -v mm="$MM" -v kin="$KINSHIP" -v idx="$INDEX" -v src="$QTL_FILE" -v max="$MAX_QTL" '
      BEGIN{r=0}
      {
        r++;
        if(max>0 && r>max) exit;

        dmr=$1; beta=$2; bsd=$3; p=$4;

        split(dmr, a, ":");
        chr_raw=a[1];
        pos=a[2];
        split(pos, b, "-"); pos=b[1];   # if format is chr:start-end

        c=chr_raw;
        sub(/^chr/,"",c);
        sub(/^ch/,"",c);
        sub(/^0+/,"",c);
        if(c=="") c=0;

        print meta, mm, kin, idx, src, r, dmr, chr_raw, c, pos, beta, bsd, p;
      }' >> "$OUT"

echo "Wrote: $OUT"
