#!/bin/bash

start_time=$(date +%s)

gzip -c 06E_subsampled_UCSC_RNA_adata.h5ad > 06E_subsampled_UCSC_RNA_adata.h5ad.gz

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

echo "Elapsed time: ${hours} hours, ${minutes} minutes, ${seconds} seconds"
