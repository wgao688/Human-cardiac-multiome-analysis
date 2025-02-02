#!/bin/bash

# Create a function to process each directory
process_dir() {
  
    # Set the genome directory and other parameters
    GENOME_DIR="/mnt/data1/william/10X_genome_refs/refdata-gex-GRCh38-2024-A/STAR_v2_7_11b_index"
    THREADS=8
    TECHNOLOGY="Dropseq"

    local dir="$1"
    # Find the R1 and R2 FASTQ files
    FASTQ_R1=$(find "$dir" -type f -name "*R1*.fastq.gz")
    FASTQ_R2=$(find "$dir" -type f -name "*R2*.fastq.gz")
    
    # Ensure both R1 and R2 are found
    if [[ -n "$FASTQ_R1" && -n "$FASTQ_R2" ]]; then
        echo "Processing directory: $dir"
        echo "R1 file: $FASTQ_R1"
        echo "R2 file: $FASTQ_R2"
        
        # Run the STARsolo script with the identified files
	bash run_STARSolo_human_heart.sh -g $GENOME_DIR -t $THREADS -1 $FASTQ_R1 -2 $FASTQ_R2 -v $TECHNOLOGY
    else
        echo "Skipping $dir: Missing R1 or R2 FASTQ files"
    fi
}

export -f process_dir  # Export the function so it can be used by xargs

# Find all immediate subdirectories and run 5 processes at a time
num_processes=5
find . -maxdepth 1 -type d | xargs -I {} -n 1 -P $num_processes bash -c 'process_dir "$@"' _ {}
