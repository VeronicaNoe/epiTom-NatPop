#!/bin/bash
# do it here: /mnt/disk2/vibanez/05_DMR-processing/05.2_DMR-annotation/ab_highPriotiryIntersection
# Define arrays
annotations=("gene-TE" "gene" "intergenic" "promotor" "TE-wo-gene")
out_file="overlap_summary.txt"

# Write header
echo -e "Annotation\tOverlap_count" > $out_file

# Loop through annotations
for annot in "${annotations[@]}"
do
    # Run intersectBed and count overlaps
    count=$(intersectBed -a ${annot}.C-DMR.methylation -b ${annot}.CG-DMR.methylation | wc -l)
    
    # Write the result
    echo -e "${annot}\t${count}" >> $out_file
done

echo "Summary written to $out_file"
