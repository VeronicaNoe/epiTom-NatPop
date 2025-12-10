#!/bin/bash

inDir="/mnt/disk6/vibanez/read_depth/bam_files"
outDir="./mean_depth_tmp"

# Create temp output directory
mkdir -p "$outDir"

# Output header
echo -e "Sample\tClean_Sample\tTissue\tMean_Depth" > mean_depth_summary_clean.tsv

# Export a function so GNU parallel can use it
process_bam () {
  bam_path="$1"
  outDir="$2"

  bam_file=$(basename "$bam_path")
  sample_name=${bam_file%.bam}

  # Clean names
  clean_sample=$(echo "$sample_name" | sed 's/_.*//')
  tissue=$(echo "$sample_name" | sed 's/^[^_]*_//' | sed 's/_Biseq\.sorted//')

  # Try to index
  if ! samtools index "$bam_path" 2>/dev/null; then
    echo "⚠️  BAM not sorted: $bam_file, re-sorting!"
    sorted_bam="${bam_path%.bam}.resorted.bam"
    samtools sort -@ 8 "$bam_path" -o "$sorted_bam"
    samtools index "$sorted_bam"
    bam_path="$sorted_bam"  # Update to re-sorted BAM
    sample_name=$(basename "$bam_path" .bam)
  fi

  # Calculate mean depth
  mean_depth=$(samtools depth "$bam_path" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')

  # Save to temp file specific to this sample
  echo -e "${sample_name}\t${clean_sample}\t${tissue}\t${mean_depth}" > "${outDir}/${sample_name}.tsv"
}

export -f process_bam

# Run all BAMs in parallel
find "$inDir" -maxdepth 1 -name "*.bam" | parallel -j 8 process_bam {} "$outDir"

# Merge all temp outputs
cat "$outDir"/*.tsv | sort >> mean_depth_summary_clean.tsv

# Clean temp directory
#rm -r "$outDir"
