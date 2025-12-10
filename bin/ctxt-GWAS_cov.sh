#!/bin/bash
start_time=$(date +%s)
##general
SAMPLE="$( cat $1 | cut -d'_' -f1,2,3 )"
COV="$( cat $1 | cut -d'_' -f4,5,6 )"
##general
#SAMPLE="$( cat $1 | cut -d'_' -f1,2 )"
#FEAT="$( cat $1 | cut -d'_' -f1)"
#COV="$( cat $1 | cut -d'_' -f3 )"
##with feature
#SAMPLE="$( cat $1 | cut -d'_' -f1,2,3 )"
#FEAT="$( cat $1 | cut -d'_' -f1,2)"
#COV="$( cat $1 | cut -d'_' -f4 )"

echo $SAMPLE $FEAT $COV
DATADIR="/mnt/disk2/vibanez/DMR-GWAS/data"
WDIR="/mnt/disk2/vibanez/descriptive-stat/ac_ctxt-GWAS"
OUTDIR="/mnt/disk2/vibanez/descriptive-stat/ac_ctxt-GWAS/results"
/mnt/disk2/vibanez/DMR-GWAS/bin/EMMAX/emmax-intel64 -v -d 10 \
	-t $DATADIR/"SNPs" \
	-p $WDIR/phenotype/$SAMPLE".phenotype" \
	-k $DATADIR/"SNPs.aBN.kinf" \
	-c $WDIR/covariables/$COV \
	-o $OUTDIR/$SAMPLE"_"$COV.cov

cat $OUTDIR/$SAMPLE"_"$COV.cov.ps | gawk '{OFS="\t"}{if($4<=5e-8){print}}' > $OUTDIR/$SAMPLE"_"$COV.cov.mQTL
pigz $OUTDIR/$SAMPLE"_"$COV.cov.ps
if [ -s $OUTDIR/$SAMPLE"_"$COV.cov.mQTL ]; then
	echo $SAMPLE "SIGNIFICATIVOS"
	mv $OUTDIR/$SAMPLE"_"$COV.cov.* $WDIR/sig
else
	echo $SAMPLE "Vacio"
        mv $OUTDIR/$SAMPLE"_"$COV.cov.* $WDIR/nonSig
fi

# End time
end_time=$(date +%s)
# Calculate the duration in seconds
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
CURRENT_DIR="$(pwd)"
rm $CURRENT_DIR/ba_targets/$SAMPLE"_"$COV.target
