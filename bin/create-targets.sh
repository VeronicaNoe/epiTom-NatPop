#!/bin/bash
cat $1 | while read line; do
  sample="$( echo $line | cut -d'|' -f1 )"

  # Find the correct directory that matches the pattern *_targets
  target_dir=$(find . -type d -name "*_targets" | head -n 1)

  if [ -n "$target_dir" ]; then
    # Create the target file and write the line to it
    touch "$target_dir/$sample.target"
    echo $line > "$target_dir/$sample.target"
  else
    echo "Target directory not found for $sample"
  fi
done
