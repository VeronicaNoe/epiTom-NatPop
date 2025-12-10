for file in *_read-depth; do
  avg_depth=$(awk '$3>0 {sum+=$3; count++} END {if (count>0) print sum/count; else print 0}' "$file");
  sample_name=$(basename "$file" | sed 's/_read-depth.*//');
  echo -e "${sample_name}\t${avg_depth}" >> read_depth_summary.tsv;
done
