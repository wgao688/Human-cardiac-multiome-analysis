#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <bam_file>"
    exit 1
fi

start_time=$(date +%s)

# Assign arguments to variables
bam_file=$1
bam_file_dir=$(dirname $bam_file)

# obtain the barcode file
#barcodes_file="$bam_file_dir/Solo.out/GeneFull/filtered/barcodes.tsv"

# Define the output directory and other parameters
output_dir=$bam_file_dir/"cellsnplite"
threads=10
min_maf=0.1
min_count=100

# Run cellsnp-lite with the specified arguments to do it at a single cell level -- this takes much longer
#cellsnp-lite -s $bam_file -b $barcodes_file -O $output_dir -p $threads --minMAF $min_maf --minCOUNT $min_count --gzip

# run cellsnp-lite at a bulk level
cellsnp-lite -s $bam_file -O $output_dir -p $threads --minMAF $min_maf --minCOUNT $min_count --cellTAG None --UMItag UB --gzip
echo "cellsnp-lite processing for $bam_file completed."

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "Elapsed time: $elapsed seconds"
