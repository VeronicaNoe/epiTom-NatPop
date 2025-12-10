#!/bin/bash
CHR="$( cat "$1" | cut -d'_' -f1)"
DMR="$( cat "$1" | cut -d'_' -f2)"
SIG="$( cat $1 | cut -d'_' -f3)"
echo $CHR $DMR $SIG
# Use all the python scripts
echo '1-Merge GWAS mQTL'
python ~/bin/merge_GWAS.py $CHR $DMR $SIG
#echo '2-Merge GWAS reml'
python ~/bin/merge_GWAS-reml.py $CHR $DMR $SIG
#echo '3-Merge GWAS mQTL from nonAssociated loci'
python ~/bin/merge_GWAS-remlNonSig.py $CHR $DMR
