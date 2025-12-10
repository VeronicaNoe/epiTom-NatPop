#!/bin/bash
DIR="/mnt/disk2/vibanez/tomato_RNAseq/read_files/"
READ="$( cat $1 | cut -d'_' -f1 )"
NAME="$( cat $1 | cut -d'_' -f2 )"

mv $DIR/$READ"_fruit_RNAseq_1.fastq.gz" $DIR/$NAME"_fruit_RNAseq_1.fastq.gz"
mv $DIR/$READ"_fruit_RNAseq_2.fastq.gz" $DIR/$NAME"_fruit_RNAseq_2.fastq.gz"
