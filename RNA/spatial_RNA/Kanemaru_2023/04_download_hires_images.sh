#!/bin/bash

## Download the hires images to each of their sample directories

input_file="03_image_metadata.txt"

# read each file, expect for the header
tail -n +2 "$input_file" | while read -r line; do

	# extract columns using awk
	sample_id=$(echo "$line" | awk '{print $1}')
	links=$(echo "$line" | awk '{print $NF}')

	output_dir="samples/$sample_id"
	for link in "${links[@]}"; do
       		wget -P "$output_dir" "$link"
		echo "$link"
	done
done
