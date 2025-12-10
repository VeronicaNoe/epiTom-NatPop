
INDIR="/mnt/disk2/vibanez/05_DMR-processing/05.1_DMR-classification/05.2_merge-DMRs/ab_pedigree/ab_merge-methylation"

# get recurrent epimutations
echo "CG-DMR recurrent epimutations"
intersectBed -a $INDIR/CG-DMR_intraG0.merged.methylation.bed \
  -b $INDIR/CG-DMR_intraG5.merged.methylation.bed -wo | wc -l

echo "C-DMR recurrent epimutations"
intersectBed -a $INDIR/C-DMR_intraG0.merged.methylation.bed \
  -b $INDIR/C-DMR_intraG5.merged.methylation.bed -wo | wc -l

