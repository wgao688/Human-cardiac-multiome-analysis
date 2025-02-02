#!/bin/bash

list_of_directories=00_directories.txt
directories=($(awk '{print $1}' $list_of_directories))

# iterate through each directory, entering it and then running the cellranger-atac count pipeline sequentially to reduce RAM usage
for directory in "${directories[@]}"; do
	echo "$directory"
	cd $directory
	cellranger-atac count --reference "/mnt/data1/william/10X_genome_refs/refdata-cellranger-arc-GRCh38-2020-A-2.0.0" --fastqs . --id counts
	cd ..
	sleep 3600
done
