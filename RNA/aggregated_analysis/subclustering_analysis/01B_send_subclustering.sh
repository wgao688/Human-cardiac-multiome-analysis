#!/bin/bash

# send the subclustering for each cell type

# run 01_run_subcluster_for_one_cell_type.py
# usage: 01_run_subcluster_for_one_cell_type.py [-h] -a ADATA_PATH -c CELL_TYPE -o OUTDIR

outdir="post_subclustering/"
mkdir -p $outdir

for adata in pre_subclustering/*.h5ad; do
	#echo $adata
	cell_type=$(basename $adata .h5ad)
	echo "$cell_type"
	nohup python3 01B_run_subcluster_for_one_cell_type.py -a $adata -c $cell_type -o $outdir &

	# wait 2 hours between jobs
	sleep 7200
done
