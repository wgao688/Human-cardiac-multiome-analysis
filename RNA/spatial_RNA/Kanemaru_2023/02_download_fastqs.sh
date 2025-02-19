#!/bin/bash

input_file="01_donor_info.txt"

# read each file, except for the header
tail -n +2 "$input_file" | while read -r line; do
	# extract columns using awk
	sample_id=$(echo "$line" | awk '{print $1}')
	fastqs=$(echo "$line" | awk '{print $NF}')

	# create new directory for each sample_id
	mkdir -p "$sample_id"

	# split the fastqs column by semicolon and download each file into the directory
	IFS=';' read -ra urls <<< "$fastqs"
	for url in "${urls[@]}"; do
		wget -P "$sample_id" "$url"
    	done
done
