#!/bin/bash
ACOL="$( cat $1| cut -d'_' -f1 )"
WIN="$( cat $1 | cut -d'_' -f2 )"
DCOL="$( cat $1 | cut -d'_' -f3 )"

OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig4/aa_identify-gene-region-affected-by-methylation/ba_get-up-down-gene-coordinates"
INDIR="/mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/aa_annotation-data"


# Extract upstream windows
echo "Processing upstream windows for $ACOL..."
awk -v win="$WIN" -v col="$ACOL" 'BEGIN {FS=OFS="\t"}
    {
        if ($4 == "+") {
            start = $2 - win;
            if (start < 1) start = 1;  # Prevent negative coordinates
            end = start + 99;
        } else if ($4 == "-") {
            end = $3 + win;
            start = end - 99;
        }
        print $1, start, end, $4, col
    }' "$INDIR/allGenes.bed" > "$OUTDIR/$ACOL.bed"

# Extract downstream windows
echo "Processing downstream windows for $DCOL..."
awk -v win="$WIN" -v col="$DCOL" 'BEGIN {FS=OFS="\t"}
    {
        if ($4 == "+") {
            start = $3 + win;
            end = start - 99;
        } else if ($4 == "-") {
            start = $2 - win;
            end = start + 99;
            if (start < 1) start = 1;  # Prevent negative coordinates
        }
        print $1, start, end, $4, col
    }' "$INDIR/allGenes.bed" > "$OUTDIR/$DCOL.bed"

echo "Upstream and downstream windows saved to $OUTDIR."
