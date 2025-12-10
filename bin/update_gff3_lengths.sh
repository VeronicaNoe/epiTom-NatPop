#!/bin/bash

# Create a mapping file with correct chromosome lengths
cat <<EOT > chrom_lengths.txt
SL2.50ch01      98543444
SL2.50ch02      55340444
SL2.50ch03      70787664
SL2.50ch04      66470942
SL2.50ch05      65875088
SL2.50ch06      49751636
SL2.50ch07      68045021
SL2.50ch08      65866657
SL2.50ch09      72482091
SL2.50ch10      65527505
SL2.50ch11      56302525
SL2.50ch12      67145203
# Add all other chromosomes and their lengths here
EOT

# Update the GFF3 file with correct sequence-region lengths
while read chrom length; do
  sed -i "s/##sequence-region ${chrom} 1 [0-9]*/##sequence-region ${chrom} 1 ${length}/" ITAG2.4_gene_models.gff3
done < chrom_lengths.txt
