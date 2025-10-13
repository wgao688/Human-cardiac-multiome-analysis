#!/bin/bash

# create a directory to store the Seurat objects
overall_Seurat_directory="Seurat_obj_dir/"
mkdir -p $overall_Seurat_directory

# iterate through each of the cell type subsetted adatas, converting them to Seurat objects
for adata in cell_type_adata/*; do
	
	#echo $adata
	cell_type=$(basename $adata .h5ad)
	echo $cell_type

	output_dir="$overall_Seurat_directory/Seurat_$cell_type/"

	# send the script to convert scanpy adata to Seurat
	nohup Rscript convert_scanpy_to_Seurat.R -a $adata -o $output_dir & 
done
