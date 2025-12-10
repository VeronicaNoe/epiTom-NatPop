# get metaplot_koDEG
GENE_LIST="/mnt/disk2/vibanez/10_data-analysis/Fig4/results"
ANNO_DIR="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"
OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig4/results"
OUTPLOT="/mnt/disk2/vibanez/10_data-analysis/Fig4/plots"
BIGWIG_DIR="/mnt/disk2/vibanez/09_KO-processing/09.3_igv-files"
# Extract coordinates for each category
for KO in ddm1 cmt3 met1; do
    echo "Processing $KO..."

grep -wFf $GENE_LIST/06.03_DEG_KO_${KO}_genes.txt $ANNO_DIR/allGenes.bed > $OUTDIR/06.04_DEG_KO-${KO}.bed
grep -wFf $GENE_LIST/06.03_nonDEG_KO_${KO}_genes.txt $ANNO_DIR/allGenes.bed > $OUTDIR/06.04_nonDEG_KO-${KO}.bed
grep -wFf $GENE_LIST/06.03_DEG_KO_${KO}_natpop_up.txt $ANNO_DIR/allGenes.bed > $OUTDIR/06.04_DEG_KO_natpop_up-${KO}.bed
grep -wFf $GENE_LIST/06.03_DEG_KO_${KO}_natpop_down.txt $ANNO_DIR/allGenes.bed > $OUTDIR/06.04_DEG_KO_natpop_down-${KO}.bed
BIGWIG_FILES=$(ls $BIGWIG_DIR/${KO}*.bw)
#
computeMatrix scale-regions \
        -S $BIGWIG_FILES \
        -R $OUTDIR/06.04_DEG_KO-$KO.bed $OUTDIR/06.04_nonDEG_KO-$KO.bed $OUTDIR/06.04_DEG_KO_natpop_up-$KO.bed $OUTDIR/06.04_DEG_KO_natpop_down-$KO.bed \
        --beforeRegionStartLength 2000 \
        --regionBodyLength 5000 \
        --afterRegionStartLength 2000 \
        --binSize 50 \
        --outFileName $OUTDIR/matrix_$KO.gz

# Generate metaplot for this KO
plotProfile -m $OUTDIR/matrix_$KO.gz \
        -out $OUTPLOT/metaplot_$KO.png \
        --perGroup \
        --numPlotsPerRow 3 \
        --colors blue red green orange purple cyan \
        --regionsLabel "DEG_KO" "nonDEG_KO" "DEG_KO & natpop up" "DEG_KO & natpop down" \
        --samplesLabel $(basename -s .bw $BIGWIG_FILES) \
        --legendLocation upper-right \
        --yMin 0 \
        --yMax auto
    echo "Finished processing $KO."
done
