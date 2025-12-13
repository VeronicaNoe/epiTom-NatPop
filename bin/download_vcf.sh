#!/bin/bash
# Download VCF from EBI

if [ ! -f SNPs_SL2.5_185-samples_biallelic_wo_indels_with_id.vcf.gz ]; then
    echo "Downloading VCF file..."
    wget -c https://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ249/ERZ24912410/SNPs_SL2.5_185-samples_biallelic_wo_indels_with_id.vcf.gz
else
    echo "File already exists. Skipping download."
fi

