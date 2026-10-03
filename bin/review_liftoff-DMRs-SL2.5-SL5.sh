#!/bin/bash
SAMPLE="$( cat $1)"
DMR_BED_GZ="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/aa_natural-accessions/ab_merge-methylation/"${SAMPLE}".merged.methylation.bed.gz"
REF_FA="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/ITAG2.4_genomic.fasta"
TGT_FA="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/SL5.0.fasta"
CHROMS="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation/chroms_SL2.5-SL5.txt"
WORK="/mnt/disk2/vibanez/otherAnalysis/01_liftoff-DMRs-SL2.5-SL5"
#Add a stable DMR ID (so we can merge back later)
zcat "$DMR_BED_GZ" | \
awk 'BEGIN{FS=OFS="\t"}
{
  id=$1 ":" $2 "-" $3;   # stable ID from SL2.5 coords
  printf "%s\t%s\t%s\t%s", $1,$2,$3,id;
  for(i=4;i<=NF;i++) printf "\t%s",$i;
  printf "\n";
}' > "$WORK/${SAMPLE}_dmrs.with_id.tsv"

#Convert 0-based BED intervals → minimal GFF3 (SL2.5 coordinates)
awk 'BEGIN{FS=OFS="\t"; print "##gff-version 3"}
{
  chr=$1; start=$2; end=$3; id=$4;

  # gene
  print chr,"DMR","gene",start,end,".",".",".","ID="id";Name="id;
  # transcript
  print chr,"DMR","mRNA",start,end,".",".",".","ID="id".t1;Parent="id;
  # exon
  print chr,"DMR","exon",start,end,".",".",".","ID="id".t1.ex1;Parent="id".t1";
}' "$WORK/${SAMPLE}_dmrs.with_id.tsv" > "$WORK/${SAMPLE}_dmrs.SL2.5.gff3"

# Run Liftoff (SL2.5 → SL5)
liftoff "$TGT_FA" "$REF_FA" \
  -g "$WORK/${SAMPLE}_dmrs.SL2.5.gff3" \
  -o "$WORK/${SAMPLE}_dmrs.SL5.gff3" \
  -u "$WORK/${SAMPLE}_dmrs.unmapped.txt" \
  -chroms "$CHROMS" \
  -p 24 \
  -a 0.90 -s 0.90

# Extract the lifted coordinates (ID → SL5 chr/start/end)
awk 'BEGIN{FS=OFS="\t"}
$3=="gene"{
  if(match($9,/ID=([^;]+)/,m)){
    id=m[1];
    print id,$1,$4,$5;   # id, chr, start, end (1-based inclusive)
  }
}' "$WORK/${SAMPLE}_dmrs.SL5.gff3" > "$WORK/${SAMPLE}_id_to_SL5.tsv"

# Rewrite your original DMR table into SL5 BED coords (0-based)
awk 'BEGIN{FS=OFS="\t"}
FNR==NR{ chr[$1]=$2; s[$1]=$3; e[$1]=$4; next }
{
  id=$4;
  if(id in chr){
    $1=chr[id]; $2=s[id]; $3=e[id];
    # drop column 4 (ID) and restore original column order
    printf "%s\t%s\t%s", $1,$2,$3;
    for(i=5;i<=NF;i++) printf "\t%s",$i;
    printf "\n";
  }
}' "$WORK/${SAMPLE}_id_to_SL5.tsv" "$WORK/${SAMPLE}_dmrs.with_id.tsv" | gzip -c > "$WORK/${SAMPLE}.methylation.SL5.coords.tsv.gz"
#
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
}' "$WORK/${SAMPLE}_id_to_SL5.tsv" "$WORK/${SAMPLE}_dmrs.with_id.tsv" | gzip -c > "$WORK/${SAMPLE}.methylation.SL5.bed.gz"
