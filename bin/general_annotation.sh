SAMPLE="$( cat $1 )"
DATA="/mnt/disk2/vibanez/02_methylkit/ai_meth-description/bb_output"
ANNOPATH="/mnt/disk2/vibanez/02_methylkit/ac_annotation/aa_data"
OUTPATH="/mnt/disk2/vibanez/02_methylkit/ai_meth-description/bd_output"
echo '------gene'
cat $DATA/$SAMPLE".meth-10mb.windows" |\
intersectBed -a - -b $ANNOPATH/"gene.anno" -wa > $OUTPATH/$SAMPLE".gene"
echo '------gene plus TE'
cat $DATA/$SAMPLE".meth-10mb.windows" |\
subtractBed -a - -b $ANNOPATH/"gene.anno" -wa |\
intersectBed -a - -b $ANNOPATH/"gene-TE.anno" -wa > $OUTPATH/$SAMPLE".gene-TE"
echo '------TE'
cat $DATA/$SAMPLE".meth-10mb.windows" | subtractBed -a - -b $ANNOPATH/"gene.anno" -wa  |\
subtractBed -a - -b $ANNOPATH/"gene-TE.anno" -wa |\
intersectBed -a - -b $ANNOPATH/"TE-wo-gene.merged.anno" -wa > $OUTPATH/$SAMPLE".TE-wo-gene"
echo '------promotor'
cat $DATA/$SAMPLE".meth-10mb.windows" | subtractBed -a - -b $ANNOPATH/"gene.anno" -wa  |\
subtractBed -a - -b $ANNOPATH/"gene-TE.anno" -wa | subtractBed -a - -b $ANNOPATH/"TE-wo-gene.merged.anno" -wa |\
intersectBed -a - -b $ANNOPATH/"promotor.merged.anno" -wa > $OUTPATH/$SAMPLE".promotor"
echo '------intergenic'
cat $DATA/$SAMPLE".meth-10mb.windows" | subtractBed -a - -b $ANNOPATH/"gene.anno" -wa |\
subtractBed -a - -b $ANNOPATH/"gene-TE.anno" -wa | subtractBed -a - -b $ANNOPATH/"TE-wo-gene.merged.anno" -wa|\
subtractBed -a - -b $ANNOPATH/"promotor.merged.anno" -wa|\
intersectBed -a - -b $ANNOPATH/"intergenic.anno" -wa > $OUTPATH/$SAMPLE".intergenic"
#echo '------exon'
cat $DATA/$SAMPLE".meth-10mb.windows" | subtractBed -a - -b $ANNOPATH/"promotor.merged.anno" -wa |\
intersectBed -a - -b $ANNOPATH/"exon-wo-TE.merged.anno" -wa > $OUTPATH/$SAMPLE".exon"
echo '------intron'
cat $DATA/$SAMPLE".meth-10mb.windows" | subtractBed -a - -b $ANNOPATH/"promotor.merged.anno" -wa | subtractBed -a - -b $ANNOPATH/"exon-wo-TE.merged.anno"-wa|\
intersectBed -a - -b $ANNOPATH/"intron-wo-TE.merged.anno" -wa > $OUTPATH/$SAMPLE".intron"

