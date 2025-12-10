#!/bin/bash

# Header
echo -e "File\tSeqs\tMin\tMax\tAvg" > fastq_stats.tsv

# Find all fastq.gz files and run seqkit stats in parallel (4 jobs at a time)
find . -maxdepth 1 -name "*.f*.gz" | xargs  -P 4 -I {} bash -c '
    FILE={}
    STATS=$(seqkit stats "$FILE" | awk "NR==2 {print \$4\"\t\"\$5\"\t\"\$6\"\t\"\$7}")
    echo -e "${FILE}\t${STATS}"
' >> fastq_stats.tsv
