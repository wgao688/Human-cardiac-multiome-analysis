#!/bin/bash

# iterate through all fragment files and tabix them
for file in "donor_fragment_files"/*.bed.gz; do
    echo "Indexing $file with tabix..."
    gunzip $file
    base_file=${file%.gz}

    # sort the bed file
    echo "Sorting $base_file..."
    sort -k1,1 -k2,2n "$base_file" -o "$base_file"

    echo "Recompressing $base_file..."
    bgzip $base_file
    echo "Indexing $base_file"
    tabix -p bed "$base_file.gz"
done

echo "all .bed.gz files indexed successfully."
