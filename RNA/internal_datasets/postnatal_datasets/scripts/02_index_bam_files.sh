#!/bin/bash

# indexes the bam files

# run the xargs parallelized version 
# -I {} tells xargs to replace {} with the filename
find . -type f -path '*/STAR/*.bam' | xargs -P 5 -I {} bash -c 'echo "{}"; samtools index "{}"'
