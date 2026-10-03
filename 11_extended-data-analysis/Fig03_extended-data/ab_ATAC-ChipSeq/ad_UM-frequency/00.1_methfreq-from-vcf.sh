#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
OUTROOT="/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq/ad_UM-frequency"
VCFDIR="/mnt/disk2/vibanez/06_get-meth-vcf/ab_epialleles"
mkdir -p "${OUTROOT}/master"
DMRSETS=(CG-DMR C-DMR)
for DMRSET in "${DMRSETS[@]}"; do
  VCF="${VCFDIR}/${DMRSET}.vcf.gz"
  [[ -s "$VCF" ]] || { echo "WARNING: missing $VCF, skipping"; continue; }
  OUT="${OUTROOT}/master/${DMRSET}.methfreq.tsv.gz"
  zcat "$VCF" \
  | awk 'BEGIN{FS=OFS="\t"}
    $0 ~ /^#/ {next}
    {
      chr=$1; pos=$2+0; end=pos+99;
      id=chr ":" pos "-" end;
      n_um=0; n_m=0; n_obs=0;
      for(i=10;i<=NF;i++){
        gt=$i;
        if(gt=="./." || gt==".|.") continue;
        n_obs++;
        if(gt=="0/0" || gt=="0|0") n_um++;
        else if(gt=="1/1" || gt=="1|1") n_m++;
      }
      fUM  = (n_obs>0 ? n_um/n_obs : "NA");
      fMeth= (n_obs>0 ? n_m/n_obs  : "NA");
      print id, chr, pos, end, n_um, n_m, n_obs, fUM, fMeth;
    }' \
  | gzip -c > "$OUT"
  echo "Wrote: $OUT"
done
