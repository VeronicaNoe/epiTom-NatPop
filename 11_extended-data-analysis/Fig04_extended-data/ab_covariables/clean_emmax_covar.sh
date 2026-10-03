#!/usr/bin/env bash
# usage: parallel -j4 'bash clean_emmax_covar.sh SV {/}' ::: SV/*.covar
#
set -euo pipefail
DIR="$1"
COVAR="$2"
MISS_FRAC_MAX="${3:-1.0}"

IN="${DIR}/${COVAR}"

[[ -s "$IN" ]] || { echo "ERROR: missing covar file: $IN" >&2; exit 1; }

BACKUP="${IN}.tmp"
REPORT="${IN}.cleanup_report.tsv"
OUT="${IN}"

cp "$IN" "$BACKUP"

awk -v OFS="\t" -v miss_frac_max="$MISS_FRAC_MAX" -v report="$REPORT" '
{
  nrow++
  nf = NF
  if(nf > max_nf) max_nf = nf

  # store whole table
  for(i=1; i<=nf; i++) cell[nrow,i] = $i

  # count rows per column
  for(i=1; i<=nf; i++) {
    total[i]++
    if($i == -9) miss[i]++
  }
}
END{
  # report header
  print "col_index", "kept", "n_rows", "n_missing", "missing_fraction" > report

  # decide columns to keep
  # always keep first 3 columns: FID IID intercept
  keep[1]=1; keep[2]=1; keep[3]=1

  for(i=4; i<=max_nf; i++){
    m = (i in miss ? miss[i] : 0)
    t = (i in total ? total[i] : nrow)
    frac = (t>0 ? m/t : 1)

    if(frac < miss_frac_max) keep[i]=1
    else keep[i]=0

    print i, keep[i], t, m, frac >> report
  }

  # also report first 3
  for(i=1; i<=3 && i<=max_nf; i++){
    m = (i in miss ? miss[i] : 0)
    t = (i in total ? total[i] : nrow)
    frac = (t>0 ? m/t : 0)
    print i, 1, t, m, frac >> report
  }

  # write cleaned file
  for(r=1; r<=nrow; r++){
    first=1
    for(i=1; i<=max_nf; i++){
      if(keep[i]){
        val = ((r,i) in cell ? cell[r,i] : "")
        if(first){ printf "%s", val; first=0 }
        else     { printf "\t%s", val }
      }
    }
    printf "\n"
  }
}
' "$BACKUP" > "$OUT"

#mv "$OUT" "$IN"

echo "Backup:  $BACKUP"
echo "Cleaned: $COVAR"
echo "Report:  $REPORT"
