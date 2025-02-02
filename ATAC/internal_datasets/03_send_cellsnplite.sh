#!/bin/bash

# For bam files in the specified directory, run cellsnp-lite on the bam files

# Function to display help message
usage() {
    echo "Usage: $0 -d <directory>"
    echo "  -d <directory>   Specify the directory to process for BAM files."
    exit 1
}

# Parse command-line options
while getopts ":d:" opt; do
    case ${opt} in
        d)
            DIR_TO_PROCESS="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

# check if the directory was provided
if [ -z "$DIR_TO_PROCESS" ]; then
    echo "Error: Directory is required." >&2
    usage
fi

cellsnplite_script="$(realpath run_cellsnplite_pseudobulked.sh)"

# change to the specified directory
cd "$DIR_TO_PROCESS" || { echo "Directory not found: $DIR_TO_PROCESS"; exit 1; }

# find and index the BAM files
n_to_process=4 # number to simultaneously process

# run this for the Fetal directories
find . -type f -path './Fetal*/counts/outs/possorted_bam.bam' | xargs -P "$n_to_process" -I {} bash -c "\"$cellsnplite_script\" \"{}\""
