#!/bin/bash
INDIR="/mnt/disk2/vibanez/GWAS/results"
OUTDIR="/mnt/disk2/vibanez/GWAS/analysis-Results/sigSNPs"
NON="/mnt/disk2/vibanez/GWAS/analysis-Results/nonSigSNPs"
SAMPLE="$( cat $1 )"
PREFIX="$( cat $1 | sed 's/.filtered.ps//g' )"
echo $SAMPLE
echo $PREFIX
if [ -s $INDIR/$PREFIX".reml" ]; then
	#zcat $INDIR/$SAMPLE | awk '{OFS="\t"}{if($4<="0.00000005"){print}}' > $OUTDIR/$PREFIX.mQTL
	cat $INDIR/$SAMPLE | awk '{OFS="\t"}{if($4<="0.00000005"){print}}' > $OUTDIR/$PREFIX".mQTL"
	if [ -s $OUTDIR/$PREFIX".mQTL" ]; then
		echo $PREFIX "SIGNIFICATIVOS"
		mv $INDIR/$PREFIX".reml"  $OUTDIR/$PREFIX".reml"
	else
		echo "Vacio"
		mv $OUTDIR/$PREFIX".mQTL" $NON/$PREFIX".mQTL"
		mv $INDIR/$PREFIX".reml"  $NON/$PREFIX".reml"
        fi
fi
