#!/usr/bin/env bash
# Usage:
#   bash 00_make_dmr_bed.sh CG-DMR_dmrs.SL5.gff3 CG-DMR_dmrs.SL5.bed
#   bash 00_make_dmr_bed.sh C-DMR_dmrs.SL5.gff3  C-DMR_dmrs.SL5.bed
DMR="$( cat $1 )"
DMR_GFF="/mnt/disk2/vibanez/otherAnalysis/01_liftoff-DMRs-SL2.5-SL5"
OUT_DIR="/mnt/disk2/vibanez/otherAnalysis/07_SV-from-graphpangenome/aa_dmrs-SL5-bed"
awk -F'\t' 'BEGIN{OFS="\t"}
  function get_attr(s, key,   n,i,a,k,v) {
    n=split(s,a,";");
    for(i=1;i<=n;i++){
      split(a[i],kv,"=");
      k=kv[1]; v=kv[2];
      if(k==key) return v;
    }
    return "NA";
  }
  $0 ~ /^#/ {next}
  $3=="gene" {
    chr=$1;
    start=$4; end=$5;

    id_sl25 = get_attr($9,"ID");                 # e.g. SL2.50ch01:13401-13500
    cov     = get_attr($9,"coverage");
    sid     = get_attr($9,"sequence_ID");
    xcn     = get_attr($9,"extra_copy_number");  # 0,1,2...
    cpid    = get_attr($9,"copy_num_ID");

    # SL5 ID from coordinates
    id_sl5 = chr ":" start "-" end;

    # BED is 0-based half-open
    bstart = start - 1;
    if(bstart < 0) bstart = 0;

    print chr, bstart, end, id_sl25, id_sl5, cov, sid, xcn, cpid;
  }' "$DMR_GFF/${DMR}_dmrs.SL5.gff3" \
| LC_ALL=C sort -k1,1V -k2,2n > "$OUT_DIR/${DMR}_dmrs.SL5.bed"
