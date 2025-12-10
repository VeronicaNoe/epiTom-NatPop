#!/bin/bash
SAMPLE="$( cat "$1" )"
OUT_FILE="${SAMPLE%.tmp}.vcf"
#OUT_DIR="/mnt/disk2/vibanez/06_get-meth-vcf/ab_epialleles"
echo $SAMPLE
Rscript --vanilla ~/bin/density_byDMR.R $SAMPLE
{
    echo "##fileformat=VCFv4.1"
    cat "$SAMPLE.tmp"
} > "$OUT_FILE"

# Remove the temporary file
rm -f "$OUT_DIR"/"$SAMPLE.tmp"
