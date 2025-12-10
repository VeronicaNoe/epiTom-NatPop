#!/bin/bash
start_time=$(date +%s)
SAMPLE="$( cat $1 | cut -d'.' -f1 )"
DMR="$( cat $1 | cut -d'_' -f2 )"
CHR="$( cat $1 | cut -d'_' -f1 )"
#echo $SAMPLE $CHR
DATADIR="/mnt/disk2/vibanez/DMR-GWAS/data"
WDIR="/mnt/disk2/vibanez/descriptive-stat/ac_ctxt-GWAS"
OUTDIR="/mnt/disk2/vibanez/descriptive-stat/ac_ctxt-GWAS/results"
/mnt/disk2/vibanez/DMR-GWAS/bin/EMMAX/emmax-intel64 -v -d 10 \
  -t $DATADIR/"SNPs" \
  -p $WDIR/phenotype/$SAMPLE".phenotype" \
  -k $DATADIR/"SNPs.aBN.kinf" \
  -o $OUTDIR/$SAMPLE

cat $OUTDIR/$SAMPLE.ps | gawk '{OFS="\t"}{if($4<=5e-8){print}}' > $OUTDIR/$SAMPLE.mQTL
pigz $OUTDIR/$SAMPLE.ps
if [ -s $OUTDIR/$SAMPLE.mQTL ]; then
	echo $SAMPLE "SIGNIFICATIVOS"
	mv $OUTDIR/$SAMPLE.* $WDIR/sig
else
	echo $SAMPLE "Vacio"
        mv $OUTDIR/$SAMPLE.* $WDIR/nonSig
fi

# End time
end_time=$(date +%s)
# Calculate the duration in seconds
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
CURRENT_DIR="$(pwd)"
rm $CURRENT_DIR/ba_targets/$SAMPLE.target
