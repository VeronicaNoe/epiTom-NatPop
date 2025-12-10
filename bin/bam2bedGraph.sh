#!/bin/bash
#set variables
SAMPLE="$( echo $1 | sed 's/.nsorted.deduplicated.bam//g' )"

context=CG
echo $SAMPLE $context
echo "track type=bedGraph" > $SAMPLE'-'${context}.bedgraph
zcat $SAMPLE.sorted.deduplicated.CX_report.txt | awk -v c=$context '$6==c && $4+$5>0 {print $1"\t"$2-1"\t"$2"\t"$4/($5+$4)}' | sort -k1,1 -k2,2n >> $SAMPLE'-'${context}.bedgraph

context=CHG
echo $SAMPLE $context
echo "track type=bedGraph" > $SAMPLE'-'${context}.bedgraph
zcat $SAMPLE.sorted.deduplicated.CX_report.txt | awk -v c=$context '$6==c && $4+$5>0 {print $1"\t"$2-1"\t"$2"\t"$4/($5+$4)}' | sort -k1,1 -k2,2n >> $SAMPLE'-'${context}.bedgraph

context=CHH
echo $SAMPLE $context
echo "track type=bedGraph" > $SAMPLE'-'${context}.bedgraph
zcat $SAMPLE.sorted.deduplicated.CX_report.txt | awk -v c=$context '$6==c && $4+$5>0 {print $1"\t"$2-1"\t"$2"\t"$4/($5+$4)}' | sort -k1,1 -k2,2n >> $SAMPLE'-'${context}.bedgraph
