#!/usr/bin/env bash
CHR="$( cat $1 )"
INDIR="/mnt/disk2/vibanez/otherAnalysis/00_public-data"
SVLEN="50"         # define "SV" as |delta| >= SVLEN (you can change later)
OUT_DIR="/mnt/disk2/vibanez/otherAnalysis/07_SV-from-graphpangenome/ab_graphpan-SV-to-bed"

zcat "$INDIR/${CHR}.vcf.gz" \
| awk -v SVLEN="$SVLEN" 'BEGIN{FS=OFS="\t"}
    /^#/ {next}
    {
      chrom=$1; pos=$2; id=$3; ref=$4; alt=$5; qual=$6; filt=$7; info=$8;

      # multi-allelic ALT: choose longest ALT for size-based classification
      n=split(alt, A, ",");
      maxAltLen=0;
      for(i=1;i<=n;i++){
        l=length(A[i]);
        if(l>maxAltLen) maxAltLen=l;
      }

      lr=length(ref);
      la=maxAltLen;

      delta = la - lr;
      absd  = (delta<0)? -delta : delta;

      # BED interval: [pos-1, pos-1+len(REF)] (>=1 bp)
      s = pos - 1;
      e = s + lr;
      if(e <= s) e = s + 1;

      cls="SNP_or_small_indel";
      if(absd >= SVLEN){
        if(delta >= SVLEN) cls="INS";
        else if(delta <= -SVLEN) cls="DEL";
        else cls="SV_other";
      }

      # Output WITHOUT sequences
      print chrom, s, e, id, pos, qual, filt, info, n, lr, la, delta, absd, cls;
    }' \
| LC_ALL=C sort -k1,1V -k2,2n \
| gzip -c > "$OUT_DIR/${CHR}.SVLEN${SVLEN}.bed.gz"


