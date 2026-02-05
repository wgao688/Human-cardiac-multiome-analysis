#!/bin/bash

# run trajectory analysis
#cell_types=(Cardiomyocyte Fibroblast Myeloid Endothelial Pericyte)
#cell_types=(Cardiomyocyte)
cell_types=(Cardiomyocyte Fibroblast Myeloid Endothelial)
#cell_types=(Pericyte)

for cell_type in ${cell_types[@]}; do
	echo "$cell_type"
	nohup python3 01_run_trajectory_per_cell_type.py -c $cell_type & 
done
