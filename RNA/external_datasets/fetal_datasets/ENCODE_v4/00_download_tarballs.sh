#!/bin/bash

# Example of a path to wget
#https://www.encodeproject.org/files/ENCFF842OKZ/@@download/ENCFF842OKZ.tar.gz

# extract the accession names and download them using the string 
accessions=($(awk '{print $5}' 00_fetal_multiome_metadata.txt | grep "ENC"))

prefix="https://www.encodeproject.org/files"
middle="@@download"
suffix=".tar.gz"

for accession in ${accessions[@]}; do
	download_path="$prefix/$accession/$middle/$accession$suffix"
	echo "$download_path"
	wget $download_path
done
