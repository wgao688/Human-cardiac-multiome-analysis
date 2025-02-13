## Create bedgraphs of specific regions 

### STEP 0: Download necessary software
# To create smaller genome windows with the bedgraph format, we need to download functions that can convert between bigwig and bedgraph
```
$ mamba install bioconda::ucsc-bedgraphtobigwig
$ mamba install bioconda::ucsc-bigwigtobedgraph
```

### STEP 1: Create bedgraph files from bigwig

The bigwig files are too big for download and for some reason run into incompatibility issues with IGV. Therefore, we will generate bedgraph files. We will specifically look at the genomic locations of some marker genes, using the same ones as for the dotplot we generated using scanpy. Run this script which will iterate through the marker genes in the file `01_marker_gene_coordinates.txt` to generate a subdirectory for each marker gene and a bedgraph file for each cell type pseudobulked. This finishes in less than 1 minute.

```
$ bash 01_create_bedgraph_files.sh 
```

As there are 11 resolved cell types, there will be 11 marker gene directories. Note that the LEC and Endocardial cells were not able to be resolved from the Endothelial cell type.

### STEP 2: Promoter only analysis

Perform the same analysis above for the promoter region only, defined as the 1000 nt upstream and downstream of the gene start. This uses `02_promoter_coordinates.txt`
```
$ bash 02_create_promoter_bedgraph_files.sh
```

### STEP 3: On laptop, generate the IGV genome browser tracks by transferring the bedgraphs over.
