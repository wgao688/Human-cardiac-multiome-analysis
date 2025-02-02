#!/bin/bash

# Set default parameter
NUM_PROCESSES=5  # Number of parallel processes

# Help function to show the arguments
show_help() {
  echo "Usage: $0 [-d DIRECTORY] [-p NUM_PROCESSES]"
  echo ""
  echo "Options:"
  echo "  -d DIRECTORY      Specify the parent directory to search (default: Fastq)"
  echo "  -s SUBDIRECTORY   Specify whether to use Gene or GeneFull directory"
  echo "  -p NUM_PROCESSES  Specify the number of parallel processes for xargs (default: 5)"
  echo "  -h                Show this help message"
}

# Parse command-line options
while getopts "d:s:p:h" opt; do
  case $opt in
    d)
      DIRECTORY=$OPTARG
      ;;
    s)
      SUBDIRECTORY=$OPTARG
      ;;
    p)
      NUM_PROCESSES=$OPTARG
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

# Ensure the directory exists
if [ ! -d "$DIRECTORY" ]; then
  echo "Error: Directory '$DIRECTORY' does not exist."
  show_help
  exit 1
fi

# Ensure that the argument is either Gene or GeneFull
if [[ $SUBDIRECTORY != "Gene" && $SUBDIRECTORY != "GeneFull" ]]; then
  echo "$SUBDIRECTORY must be either Gene (only the exonic counts) or GeneFull (intronic + exonic counts)" 
fi	

# Find directories that contain count matrices in "Solo.out/GeneFull" which corresponds
# specify maxdepth = 4 to traverse down 4 directories in $DIRECTORY/<sample>/STAR/GeneFull
count_matrix_paths=$(find "$DIRECTORY" -maxdepth 4 -type d -name $SUBDIRECTORY -print)

# Check if count_matrix_paths is not empty
if [ -z "$count_matrix_paths" ]; then
  echo "No count matrix paths found in the specified directory."
  exit 0
fi

# Process each path using xargs to run multiple processes in parallel
DIR=$(dirname $0)
run_SoupX_script=$DIR/run_SoupX.R
echo "$count_matrix_paths" | xargs -I {} -n 1 -P $NUM_PROCESSES bash -c '
  path="$1"
  echo "Processing count matrix path: $path"
  Rscript '$run_SoupX_script' "$path"
' _ {}

echo "All tasks submitted."
