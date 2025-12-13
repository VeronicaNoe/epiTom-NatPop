#!/bin/bash
echo $1
SNP_NAME="$( cat $1 | cut -d'|' -f1 )"
SAMPLE="$( cat $1 | cut -d'|' -f2 )"

REF_FASTA="02_pseudo-pan-genomes/02.0_SL2.50-genome/ITAG2.4_genomic.fasta"
REF_SNP="01_raw_data/01.3_vcfiles/SNPs_*.gz"
OUT_DIR="02_pseudo-pan-genomes/02.1_make-pseudogenomes"
echo "---------SAMPLE: " $SAMPLE
mkdir -p ${OUT_DIR}/$SAMPLE
bcftools consensus ${REF_SNP} --sample $SNP_NAME --fasta-ref ${REF_FASTA} > ${OTU_DIR}/$SAMPLE/$SAMPLE"_pseudoPangenome.fa"
