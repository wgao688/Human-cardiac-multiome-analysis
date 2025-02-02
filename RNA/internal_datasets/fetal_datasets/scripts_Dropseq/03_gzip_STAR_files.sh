#!/bin/bash

# gzip the STAR output files, since this is the format the SoupX and scanpy expect

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

# Check if the directory was provided
if [ -z "$DIR_TO_PROCESS" ]; then
    echo "Error: Directory is required." >&2
    usage
fi

# Change to the specified directory
cd "$DIR_TO_PROCESS" || { echo "Directory not found: $DIR_TO_PROCESS"; exit 1; }

nohup gzip $DIR_TO_PROCESS/*/STAR/Solo.out/Gene/filtered/* &
nohup gzip $DIR_TO_PROCESS/*/STAR/Solo.out/GeneFull/filtered/* &
nohup gzip $DIR_TO_PROCESS/*/STAR/Solo.out/Gene/raw/* &
nohup gzip $DIR_TO_PROCESS/*/STAR/Solo.out/GeneFull/raw/* &
