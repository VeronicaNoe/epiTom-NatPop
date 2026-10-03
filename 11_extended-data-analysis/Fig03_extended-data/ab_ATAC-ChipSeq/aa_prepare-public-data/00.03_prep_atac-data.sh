OUTROOT="/mnt/disk2/vibanez/otherAnalysis/05_ATAC-ChipSeq/aa_prepare-public-data"
BW="${OUTROOT}/atlas_bw"
OUTTAB="${OUTROOT}/master/CG-DMR.ATAC.mean_by_DMR.tsv"

BED4_GZ="/mnt/disk2/vibanez/otherAnalysis/04_liftoff-DMRs-SL2.5-M82/CG-DMR.M82.hiconf.mapq30.chr.bed4.gz"
BED4="${OUTROOT}/master/CG-DMR.M82.bed4"
mkdir -p "${OUTROOT}/master"
zcat "$BED4_GZ" > "$BED4"

ATAC1=$(ls -1 "${BW}"/*ATAC*M82*rep1*.bigwig 2>/dev/null | head -n1)
ATAC2=$(ls -1 "${BW}"/*ATAC*M82*rep2*.bigwig 2>/dev/null | head -n1)

multiBigwigSummary BED-file \
  --bwfiles "$ATAC1" "$ATAC2" \
  --BED "$BED4" \
  --labels ATAC_rep1 ATAC_rep2 \
  --outFileName "${OUTROOT}/master/CG-DMR.ATAC.mean_by_DMR.npz" \
  --outRawCounts "$OUTTAB"
