#!/bin/bash

start_time=$(date +%s)

gzip -c 06C_final_adata_for_UCSC.h5ad > 06C_final_adata_for_UCSC.h5ad.gz

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

echo "Elapsed time: ${hours} hours, ${minutes} minutes, ${seconds} seconds"
