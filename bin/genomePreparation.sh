#!/bin/bash
SAMPLE="$( cat "$1" )"
echo $SAMPLE
WD="/mnt/disk2/vibanez/03_biseq-processing/03.0_genome-preparation"
bismark_genome_preparation --bowtie2 --verbose $WD/$SAMPLE
