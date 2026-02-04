#!/bin/bash

# sender script to run silhouette analysis for multiple cell types

# list of cell types
#cell_types=(
#  "Adipocyte" "Cardiomyocyte" "Endocardial" "Endothelial" "Epicardial"
#  "Fibroblast" "Lymphoid" "LEC" "Mast" "Myeloid" "Neuronal" "Pericyte" "vSMC"
#)

cell_types=("Pericyte")

# directory for output logs
outdir="logs_silhouette/"
mkdir -p $outdir

# iterate through cell types
for cell_type in "${cell_types[@]}"; do
    echo "Submitting silhouette analysis for $cell_type"

    nohup python3 05_run_silhouette_single_cell_type.py -c $cell_type \
        > "${outdir}/${cell_type}.out" 2> "${outdir}/${cell_type}.err" &
done

echo "All jobs submitted!"
