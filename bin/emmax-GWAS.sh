#!/bin/bash
start_time=$(date +%s)
SAMPLE="$( cat $1 | cut -d'.' -f1 )"
DMR="$( cat $1 | cut -d'_' -f2 )"
CHR="$( cat $1 | cut -d'_' -f1 )"
echo $SAMPLE $CHR
WDIR="/mnt/disk2/vibanez/GWAS"
OUTDIR="/mnt/disk2/vibanez/GWAS/results"
/home/vibanez/bin/EMMAX/emmax-intel64 -v -d 10 \
  -t $WDIR"/data/SNPs_SL2-5.LD" \
  -p $WDIR"/phenotype/"$SAMPLE".phenotype" \
  -k $WDIR"/data/SNPs_SL2-5.LD.aBN.kinf" \
  -o $OUTDIR/$SAMPLE

cat $OUTDIR/$SAMPLE.ps | awk '{OFS="\t"}{if($4<=5e-8){print}}' > $OUTDIR/$SAMPLE.mQTL
pigz $OUTDIR/$SAMPLE.ps
mv  $OUTDIR/$SAMPLE.ps.gz $WDIR/toBrick
if [ -s $OUTDIR/$SAMPLE.mQTL ]; then
	echo $SAMPLE "SIGNIFICATIVOS"
	mv $OUTDIR/$SAMPLE.* $WDIR/$CHR/sig
else
	echo $SAMPLE "Vacio"
        mv $OUTDIR/$SAMPLE.* $WDIR/$CHR/NonSig
fi

# End time
end_time=$(date +%s)
# Calculate the duration in seconds
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
rm $WDIR/ba_targets/$SAMPLE.mQTL.target
