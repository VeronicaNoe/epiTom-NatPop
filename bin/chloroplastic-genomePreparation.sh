#!/bin/bash
INDIR="/mnt/disk2/vibanez/03_biseq-processing/03.6_conversion-rate/aa_genome-preparation"
bismark_genome_preparation --bowtie2 --verbose ${INDIR}
