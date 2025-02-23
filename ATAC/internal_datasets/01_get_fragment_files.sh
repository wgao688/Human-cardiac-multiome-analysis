#!/bin/bash

# create fragment_files directory
mkdir -p fragment_files

# iterate through all subdirectores in raw_data/
for dir in raw_data/; do
	# copy the fragments files to the fragment_files/ directory using the subdirectory name
	cp "$dir/counts/outs/fragments.tsv.gz" "fragment_files/${dir%/}_fragments.tsv.gz"
	# also copy the index files
	cp "$dir/counts/outs/fragments.tsv.gz.tbi" "fragment_files/${dir%/}_fragments.tsv.gz.tbi"
	echo "Copied fragments for $dir"
done
