INDIR="/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/bd_get-DMR-meth-gene-expression-correlation/tmp"
OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/be_merged_correlations"
OUTPUT_FILE="$OUTDIR/allCorrelation_meth-geneExpression.tsv"
ANNO_DIR=""
# Combine files while skipping headers
find "$INDIR" -maxdepth 1 -name "*correlation" -exec awk 'FNR > 1' {} + > "$OUTPUT_FILE"
cat allCorrelation_meth-geneExpression.tsv |\
awk '{OFS="\t"}{print $1, $2, $2+99,$3, $4, $6}' |\
intersectBed -a - -b ${ANNO_DIR}/TE-wo-gene.anno -wa > DMRs_over_TEs.bed
