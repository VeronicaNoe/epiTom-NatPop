WIN="$( cat $1 )"
WINPATH="/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/ba_get-up-down-gene-coordinates"
OUTPATH="/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/bb_get-meth-levels-around-gene-windows"
DMRPATH="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ae_mergedAnnotation"

intersectBed -a $WINPATH/$ANNO.bed -b $DMRPATH/general.C-DMR.methylation -wb > $OUTPATH/$WIN.C-DMR
intersectBed -a $WINPATH/$ANNO.bed -b $DMRPATH/general.CG-DMR.methylation -wb > $OUTPATH/$WIN.CG-DMR
