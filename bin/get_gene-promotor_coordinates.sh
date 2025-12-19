WIN="$( cat $1 )"
WINPATH="10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/ba_get-up-down-gene-coordinates"
OUTPATH="10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/bb_get-meth-levels-around-gene-windows"
DMRPATH="05_DMR-processing/05.2_DMR-annotation/ae_mergedAnnotation"

intersectBed -a $WINPATH/$ANNO.bed -b $DMRPATH/general.C-DMR.methylation -wb > $OUTPATH/$WIN.C-DMR
intersectBed -a $WINPATH/$ANNO.bed -b $DMRPATH/general.CG-DMR.methylation -wb > $OUTPATH/$WIN.CG-DMR
