#!/bin/bash

# indexes the bam files in the subdirectories of the given directory

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

# change to the specified directory
cd "$DIR_TO_PROCESS" || { echo "Directory not found: $DIR_TO_PROCESS"; exit 1; }

# find and index the BAM files
find . -type f -path '*/STAR/*.bam' | xargs -P 5 -I {} bash -c 'echo "{}"; samtools index "{}"'
