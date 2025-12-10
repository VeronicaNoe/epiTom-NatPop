#!/bin/bash
SAMPLE="$( echo $1 | cut -d'_' -f1,2)"
INDIR="/mnt/disk2/vibanez/02_methylkit/am_KO-RNAseq/sorted-bam"
echo "========  SAMPLE  ========="
echo $SAMPLE
echo $1
echo "========= bam 2 bg ========"
bedtools genomecov -ibam $INDIR/$SAMPLE"_RNAseq.sorted.bam" -bg | grep 'SL2.50ch' > $INDIR/$SAMPLE"_RNAseq.bg"
#"========= bg 2 bw   ========="
#bedGraphToBigWig e.bedGraph chrom.sizes yourfile.bw
