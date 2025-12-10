#!/bin/bash
input_file=$1

# run R.r script with ${_output_file} parameter
Rscript --vanilla ../../00_DMR_methylkit_for_parallel.R ${input_file} 

echo
echo "----- DONE -----"
echo

exit 0
