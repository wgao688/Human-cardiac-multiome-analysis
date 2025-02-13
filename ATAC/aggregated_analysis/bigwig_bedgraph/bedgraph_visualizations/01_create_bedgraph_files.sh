#!/bin/bash

# for each of the marker genes, create bedgraphs per cell type for each marker gene
# create a directory for each marker gene

# save the marker genes and their genomic locations in arrays
marker_gene_file="01_marker_gene_coordinates.txt"
gene_names=($(grep -v "#" $marker_gene_file | awk '{print $1}'))
num_genes=${#gene_names[@]}
cell_types=($(grep -v "#" $marker_gene_file | awk '{print $2}'))

chromosome_locs=($(grep -v "#" $marker_gene_file | awk '{print $3}'))
start_locs=($(grep -v "#" $marker_gene_file | awk '{print $4}'))
end_locs=($(grep -v "#" $marker_gene_file | awk '{print $5}'))

echo "number of marker genes: $num_genes"

# iterate through all of the marker genes
for i in $(seq 0 $((num_genes-1)) ); do
	gene_name=${gene_names[i]}
	chr=${chromosome_locs[i]}
	start_pos=${start_locs[i]}
	end_pos=${end_locs[i]}

	echo "Creating bedfiles for $gene_name: $chr:$start_pos-$end_pos"

	# create directory
	mkdir -p whole_gene/$gene_name

	# generate a bigwig file for each cell type
	for j in $(seq 0 $((num_genes-1)) ); do
		cell_type=${cell_types[j]}
		echo "$cell_type"
		bigWigToBedGraph -chrom=$chr -start=$start_pos -end=$end_pos ../cell_type_bigwig/$cell_type.bw whole_gene/"$gene_name"/"$cell_type"_"$gene_name".bedGraph
	done
done
