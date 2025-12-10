#!/bin/bash
cat $1 | awk '{ print $1, $2, $3, $6, $7, $8 }' | sed 's/ /\t/g' > ../../af_merge_C-DMRs/aa_data/$1
