#!/bin/bash

# Function to display help message
show_help() {
  echo "Usage: $0 [-d DIRECTORY] [-h]"
  echo ""
  echo "Options:"
  echo "  -d DIRECTORY    Specify the directory to search for .fastq.gz files"
  echo "  -h              Show help message"
}

# Parse options
while getopts "d:h" opt; do
  case $opt in
    d)
      directory=$OPTARG
      ;;
    h)
      show_help
      exit 0
      ;;
    *)
      show_help
      exit 1
      ;;
  esac
done


# Loop over each .fastq.gz file in the specified directory

# Check if directory is provided, otherwise exit
if [[ -z "$directory" ]]; then
  echo "Error: No directory specified."
  show_help
  exit 1
fi

# Change to the specified directory
if [[ ! -d "$directory" ]]; then
  echo "Error: The directory '$directory' does not exist."
  exit 1
fi

cd "$directory" || exit

for file in *.fastq.gz; do
    # Extract the directory name from the filename before "_S*_" using pattern matching, as this corresponds to the sample name
    dir_name=$(echo "$file" | awk -F'_S[0-9]+' '{print $1}')

    echo "$dir_name"
    # Create the directory if it doesn't already exist
    mkdir -p "$dir_name"

    # Move the file into the corresponding directory
    mv "$file" "$dir_name/"
done

echo "Files have been organized into directories, one per sample"
