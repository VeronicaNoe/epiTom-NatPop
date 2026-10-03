#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

DMRSET="$( cat $1 )"
MAPQ_MIN="${MAPQ_MIN:-30}"   # override: MAPQ_MIN=20 bash ...

WORK="/mnt/disk2/vibanez/otherAnalysis/01_liftoff-DMRs-SL2.5-M82"
cd "$WORK"

INTER="${WORK}/${DMRSET}_intermediate"
HICONF="${WORK}/${DMRSET}_id_to_M82.hiconf.tsv"
WITHID="${WORK}/${DMRSET}_dmrs.with_id.tsv"

OUT_MAPQ="${WORK}/${DMRSET}_id_to_mapq.clean.tsv"
OUT_HI_MAPQ="${WORK}/${DMRSET}_id_to_M82.hiconf.mapq${MAPQ_MIN}.tsv"
OUT_BED4="${WORK}/${DMRSET}.M82.hiconf.mapq${MAPQ_MIN}.bed4.gz"
OUT_BED4_CHR="${WORK}/${DMRSET}.M82.hiconf.mapq${MAPQ_MIN}.chr.bed4.gz"

# 1) Compute id -> MAPQ from SAMs (primary only)
gawk 'BEGIN{FS=OFS="\t"}
$0!~/^@/{
  # extract base DMR id like SL2.50ch01:13401-13500
  if(match($1,/SL2\.50ch[0-9]+:[0-9]+-[0-9]+/,m)) q=m[0]; else next;
  flag=$2+0;
  if(and(flag,256)==0 && and(flag,2048)==0){
    if($5 > mq[q]) mq[q]=$5;
  }
}
END{for(q in mq) print q,mq[q]}' \
"${INTER}"/SL2.50ch*_to_chr*.sam > "$OUT_MAPQ"

# 2) Filter hiconf by MAPQ threshold
awk -v MQ="$MAPQ_MIN" 'BEGIN{FS=OFS="\t"}
FNR==NR{mq[$1]=$2; next}
{
  id=$1;
  if((id in mq) && mq[id] >= MQ) print $1,$2,$3,$4;
}' "$OUT_MAPQ" "$HICONF" > "$OUT_HI_MAPQ"

# 3) Make M82 BED4 with id (GFF 1-based inclusive -> BED 0-based end-exclusive)
awk 'BEGIN{FS=OFS="\t"}
FNR==NR{chr[$1]=$2; s[$1]=$3; e[$1]=$4; next}
{
  id=$4;
  if(id in chr){
    bed_start = s[id]-1;
    bed_end   = e[id];
    print chr[id], bed_start, bed_end, id;
  }
}' "$OUT_HI_MAPQ" "$WITHID" | gzip -c > "$OUT_BED4"

# 4) Restrict to canonical chr1-12 (keeps you consistent with mark files)
zcat "$OUT_BED4" | awk '$1 ~ /^chr([1-9]|1[0-2])$/' | gzip -c > "$OUT_BED4_CHR"

# 5) QC
echo "== ${DMRSET} =="
echo "hiconf:"        $(wc -l < "$HICONF")
echo "hiconf_mapq:"   $(wc -l < "$OUT_HI_MAPQ")
echo "bed4_chr:"      $(zcat "$OUT_BED4_CHR" | wc -l)

echo "length dist (chr1-12 BED4):"
zcat "$OUT_BED4_CHR" | awk '{print $3-$2}' | sort -n | uniq -c | head
