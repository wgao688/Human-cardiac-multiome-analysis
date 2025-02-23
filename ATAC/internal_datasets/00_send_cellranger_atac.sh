#!/bin/bash

raw_data_dir="raw_data/"
list_of_directories=00_directories.txt
directories=($(awk '{print $1}' list_of_directories))
reference_genome_path="/mnt/data1/william/10X_genome_refs/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"


# iterate through each directory, entering it and then running the cellranger-atac count pipeline serially to reduce RAM usage
# this will produce a fragment file in $DIR/counts/outs/

for directory in "${directories[@]}"; do
	
	echo "$directory"
	
	cd $directory

	cellranger-atac count --reference $reference_genome_path --fastqs . --id counts
	cd ..

	# wait for 1 hour
	sleep 3600
done
