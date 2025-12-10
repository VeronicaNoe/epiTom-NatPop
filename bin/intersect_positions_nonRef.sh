INDIR="/mnt/disk2/vibanez/02_methylkit/nonRefSeq/ab_merged"
ANNODIR="/mnt/disk2/vibanez/02_methylkit/nonRefSeq/ac_annotation/aa_data"
OUTDIR="/mnt/disk2/vibanez/02_methylkit/nonRefSeq/ac_annotation/positions"

echo '------gene'
zcat $INDIR/CG_methylation.filtered.gz | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNODIR/"nonRef-genes.anno"  > $OUTDIR/CG_gene.bed
zcat $INDIR/CG_methylation.filtered.gz | awk '{OFS="\t"}{if (NR>1){print $1, $2, $3} }' |\
intersectBed -a - -b $ANNODIR/"nonRef-genes.anno" -wao | awk '{if($7>0){print $4,$5,$6}}' > $OUTDIR/CG_gene.ID

zcat $INDIR/CHG_methylation.filtered.gz | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNODIR/"nonRef-genes.anno"  > $OUTDIR/CHG_gene.bed
zcat $INDIR/CHG_methylation.filtered.gz | awk '{OFS="\t"}{if (NR>1){print $1, $2, $3} }' |\
intersectBed -a - -b $ANNODIR/"nonRef-genes.anno" -wao | awk '{if($7>0){print $4,$5,$6}}' > $OUTDIR/CHG_gene.ID

zcat $INDIR/CHH_methylation.filtered.gz | awk '{OFS="\t"}{if (NR>1){print} }' |\
intersectBed -a - -b $ANNODIR/"nonRef-genes.anno"  > $OUTDIR/CHH_gene.bed
zcat $INDIR/CHH_methylation.filtered.gz | awk '{OFS="\t"}{if (NR>1){print $1, $2, $3} }' |\
intersectBed -a - -b $ANNODIR/"nonRef-genes.anno" -wao | awk '{if($7>0){print $4,$5,$6}}' > $OUTDIR/CHH_gene.ID
