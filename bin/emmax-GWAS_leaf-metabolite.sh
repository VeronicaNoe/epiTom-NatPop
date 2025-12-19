#!/bin/bash
start_time=$(date +%s)
META="$( cat $1 | cut -d'_' -f1 )"
META_DIR="10_data-analysis/Fig5/aa_GWAS-metabolome/bb_phenotype"

OUT_DIR="10_data-analysis/Fig5/aa_GWAS-metabolome/bd_results"
MM_DIR="10_data-analysis/Fig5/aa_GWAS-metabolome/ba_markers"
MM="$( cat $1 | cut -d'_' -f2 )"

KINSHIP="$( cat $1 | cut -d'_' -f3 )"
INDEX="$( cat $1 | cut -d'_' -f4 )"
OUTPUT="${META}_${MM}_${KINSHIP}_${INDEX}"
# Count the number of lines in the .tped file
NUM_LINES=$(wc -l < "${MM_DIR}/${MM}_general_leaf_LD.tped" )
#THRESHOLD=$(awk -v n="$NUM_LINES" 'BEGIN {print 0.05 / n}')
if [[ "$MM" == "DMR" ]]; then
    THRESHOLD="1.048466e-07"
else
    THRESHOLD="0.000000041"
fi
echo "${THRESHOLD}"

~/bin/tools/EMMAX/emmax-intel64 -v -d 10 \
  -t ${MM_DIR}/${MM}"_general_leaf_LD" \
  -p ${META_DIR}/${META} \
  -k ${MM_DIR}/${KINSHIP}"_general_leaf_LD."${INDEX}".kinf" \
  -o ${OUT_DIR}/${OUTPUT}

cat ${OUT_DIR}/${OUTPUT}.ps |\
awk -v threshold="$THRESHOLD" '{OFS="\t"} {if ($4 <= threshold) {print}}' > ${OUT_DIR}/${OUTPUT}.QTL
pigz ${OUT_DIR}/${OUTPUT}.ps

if [ -s ${OUT_DIR}/${OUTPUT}.QTL ]; then
        echo $SAMPLE "SIGNIFICATIVOS"
        mv ${OUT_DIR}/${OUTPUT}.* ${OUT_DIR}/sig/
else
        echo $SAMPLE "Non-sig"
        mv ${OUT_DIR}/${OUTPUT}.* ${OUT_DIR}/nonSig/
fi

# End time
end_time=$(date +%s)
# Calculate the duration in seconds
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
