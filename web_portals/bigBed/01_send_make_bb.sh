#!/bin/bash

ATAC_fragment_dir="/mnt/data1/william/human_heart_project/Final_manuscript_analysis/web_portals/cellxgene/ATAC_fragment_files"

for bed in $ATAC_fragment_dir/*.bed.gz; do
	# send the script to create a bigBed file from the fragment files
	nohup bash 01_make_bigbed.sh -b $bed & 
	sleep 5
done
