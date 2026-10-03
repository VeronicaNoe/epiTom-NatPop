#!/usr/bin/env bash
CHR="$(cat "$1")"
INDIR="/mnt/disk2/vibanez/otherAnalysis/10_unique-SV-from-graphpangenome/aa_graphpan-unique-vcf"
OUT_DIR="/mnt/disk2/vibanez/otherAnalysis/10_unique-SV-from-graphpangenome/ab_unique-graphpan-SV-to-bed"
GENOME="/mnt/disk2/vibanez/otherAnalysis/07_SV-from-graphpangenome/SL5.chrom.sizes"

SVLEN="${2:-50}"
mkdir -p "$OUT_DIR"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' "$INDIR/${CHR}.vcf.gz" \
| awk -v SVLEN="$SVLEN" 'BEGIN{FS=OFS="\t"}
    {
      chrom=$1; pos=$2; id=$3; ref=$4; alt=$5; qual=$6; filt=$7; info=$8;

      nALT = split(alt, A, ",");
      maxAltLen=0;
      for(i=1;i<=nALT;i++){
        l=length(A[i]);
        if(l>maxAltLen) maxAltLen=l;
      }

      lr=length(ref);
      la=maxAltLen;

      delta = la - lr;
      absd  = (delta<0)? -delta : delta;

      s = pos - 1;
      e = s + lr;
      if(e <= s) e = s + 1;

      cls="SNP_or_small_indel";
      if(absd >= SVLEN){
        if(delta > 0) cls="INS";
        else if(delta < 0) cls="DEL";
        else cls="SV_other";
      }

      # chr start end id pos qual filter info nALT lenREF lenALT delta absdelta vclass
      print chrom, s, e, id, pos, qual, filt, info, nALT, lr, la, delta, absd, cls;
    }' \
| bedtools sort -g "$GENOME" -i - \
| bgzip -c > "$OUT_DIR/${CHR}.SVLEN${SVLEN}.bed.gz"

tabix -p bed "$OUT_DIR/${CHR}.SVLEN${SVLEN}.bed.gz"
