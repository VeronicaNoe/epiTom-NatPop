#!/bin/bash
SAMPLE="$( cat $1 )"
echo $SAMPLE
#Rscript --vanilla ~/bin/split_SNPs-per-DMRs.R $SAMPLE
#python ~/bin/process_chunks.py $SAMPLE  --chunk_size 10000 --num_cpus 45
python ~/bin/combine_GWAS-pvalues-chisq.py $SAMPLE
#python ~/bin/combine_GWAS-pvalues.py $SAMPLE
#python ~/bin/combine_GWAS-chisq.py $SAMPLE
