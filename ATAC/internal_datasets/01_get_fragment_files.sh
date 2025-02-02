#!/bin/bash

# create fragment_files directory if it doesn't exist
mkdir -p fragment_files

# Loop through all subdirectories
for dir in */; do
  # Check if the subdirectory contains 'counts/outs/' directory
  if [ -d "$dir/counts/outs/" ]; then
    # Check if the fragments.tsv.gz and fragments.tsv.gz.tbi files exist
    if [ -f "$dir/counts/outs/fragments.tsv.gz" ] && [ -f "$dir/counts/outs/fragments.tsv.gz.tbi" ]; then
      # Copy the fragments files to the fragment_files/ directory using the subdirectory name
      cp "$dir/counts/outs/fragments.tsv.gz" "fragment_files/${dir%/}_fragments.tsv.gz"
      cp "$dir/counts/outs/fragments.tsv.gz.tbi" "fragment_files/${dir%/}_fragments.tsv.gz.tbi"
      echo "Copied fragments for $dir"
    else
      echo "Fragments files not found in $dir"
    fi
  else
    echo "'counts/outs/' not found in $dir"
  fi
done
