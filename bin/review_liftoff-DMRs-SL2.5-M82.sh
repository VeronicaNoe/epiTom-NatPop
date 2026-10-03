#!/bin/bash
set -euo pipefail
SAMPLE="$(cat "$1")"
DMR_BED_GZ="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation/${SAMPLE}.merged.methylation.bed.gz"
REF_FA="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/ITAG2.4_genomic.fasta"
TGT_FA="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/SLM82.fasta"
CHROMS="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/chroms_SL2.5-M82.txt"
WORK="/mnt/disk2/vibanez/otherAnalysis/04_liftoff-DMRs-SL2.5-M82"

# 1) Add stable DMR ID (SL2.5 coords)
zcat "$DMR_BED_GZ" | awk 'BEGIN{FS=OFS="\t"}
{
  id=$1 ":" $2 "-" $3;
  printf "%s\t%s\t%s\t%s", $1,$2,$3,id;
  for(i=4;i<=NF;i++) printf "\t%s",$i;
  printf "\n";
}' > "$WORK/${SAMPLE}_dmrs.with_id.tsv"

# 2) BED-like intervals -> minimal GFF3
# IMPORTANT:
# - If your DMR file is 1-based inclusive (likely), keep start=$2; end=$3 (as below)
# - If your DMR file is true BED (0-based, end-exclusive), change to: start=$2+1; end=$3
awk 'BEGIN{FS=OFS="\t"; print "##gff-version 3"}
{
  chr=$1; start=$2; end=$3; id=$4;

  print chr,"DMR","gene",start,end,".",".",".","ID="id";Name="id;
  print chr,"DMR","mRNA",start,end,".",".",".","ID="id".t1;Parent="id;
  print chr,"DMR","exon",start,end,".",".",".","ID="id".t1.ex1;Parent="id".t1";
}' "$WORK/${SAMPLE}_dmrs.with_id.tsv" > "$WORK/${SAMPLE}_dmrs.SL2.5.gff3"

# 3) Liftoff SL2.5 -> M82
liftoff "$TGT_FA" "$REF_FA" \
  -g "$WORK/${SAMPLE}_dmrs.SL2.5.gff3" \
  -o "$WORK/${SAMPLE}_dmrs.M82.gff3" \
  -u "$WORK/${SAMPLE}_dmrs.unmapped.txt" \
  -chroms "$CHROMS" \
  -exclude_partial \
  -dir "$WORK/${SAMPLE}_intermediate" \
  -p 16 \
  -a 0.90 -s 0.90

# 4) Extract lifted coords + (sequence_ID, coverage) from gene lines
awk 'BEGIN{FS=OFS="\t"}
$3=="gene"{
  id=""; seqid="NA"; cov="NA";

  if(match($9,/ID=([^;]+)/,m)) id=m[1];
  if(match($9,/sequence_ID=([^;]+)/,m)) seqid=m[1];
  if(match($9,/coverage=([^;]+)/,m)) cov=m[1];

  if(id!="") print id,$1,$4,$5,seqid,cov;   # id, chr, gff_start, gff_end, seqID, cov
}' "$WORK/${SAMPLE}_dmrs.M82.gff3" > "$WORK/${SAMPLE}_id_to_M82.full.tsv"

# 5) High-confidence filter (adjust thresholds as you like)
# Here: keep cov>=0.95 and sequence_ID>=0.95
awk 'BEGIN{FS=OFS="\t"}
{
  id=$1; chr=$2; s=$3; e=$4; seqid=$5+0; cov=$6+0;
  if(seqid>=0.95 && cov>=0.95) print id,chr,s,e;
}' "$WORK/${SAMPLE}_id_to_M82.full.tsv" > "$WORK/${SAMPLE}_id_to_M82.hiconf.tsv"

# 6) Write M82 BED (0-based) for downstream bedtools intersects
awk 'BEGIN{FS=OFS="\t"}
FNR==NR{ chr[$1]=$2; s[$1]=$3; e[$1]=$4; next }
{
  id=$4;
  if(id in chr){
    bed_start = s[id]-1;
    bed_end   = e[id];
    printf "%s\t%d\t%d", chr[id], bed_start, bed_end;
    for(i=5;i<=NF;i++) printf "\t%s",$i;
    printf "\n";
  }
}' "$WORK/${SAMPLE}_id_to_M82.hiconf.tsv" "$WORK/${SAMPLE}_dmrs.with_id.tsv" \
| gzip -c > "$WORK/${SAMPLE}.methylation.M82.bed.gz"

echo "TOTAL:"   $(wc -l < "$WORK/${SAMPLE}_dmrs.with_id.tsv")
echo "MAPPED_hiConf:" $(wc -l < "$WORK/${SAMPLE}_id_to_M82.hiconf.tsv")
echo "UNMAPPED:" $(wc -l < "$WORK/${SAMPLE}_dmrs.unmapped.txt")
