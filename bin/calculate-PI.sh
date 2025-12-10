#!/bin/bash
GROUP="$( cat $1 | cut -d'_' -f1 )"
MARKER="$( cat $1 | cut -d '_' -f2 )"
ANNOTATION="$( cat $1 | cut -d'_' -f3 )"
SAMPLE="$( cat $1 )"
echo $SAMPLE
OUTDIR="/mnt/disk2/vibanez/10_data-analysis/Fig2/ac_calculate-nucleotide-dmr-diversity/ba_get-PI"
INDIR="/mnt/disk2/vibanez/01_raw-data/01.2_vcfiles/ac_processed-vcf/ba_split-by-groups"
#INDIR="/mnt/data6/vibanez/SNPs/vcfiles/aa_split-vcf-group/bb_output"
# check we did all of this with 20 samples because we need to get balance number of acc per group
vcftools --gzvcf $INDIR/${SAMPLE}".vcf.gz" --thin 1 --window-pi 100 --out $OUTDIR/$SAMPLE
#parallel -j70 'vcftools --vcf {} --thin 0 --window-pi 100 --out /{/.}' ::: .vcf}
