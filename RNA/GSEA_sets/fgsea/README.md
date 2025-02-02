### Reducing redundant gene sets via Jaccard similarity

This directory goes over how reduce the redundancy of gene sets by clustering based on Jaccard similarity.

### STEP 0: Download the gene sets from here:

To download the Gene Ontology Biological Processes (GO:BP) set
```
wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/c5.go.bp.v2023.2.Hs.symbols.gmt 
```

### STEP 1: After this, remove redundant sets by computing the Jaccard distance using this script between all GO terms. This represents a degree of overlap and stores the results as a distance matrix, and then performs hierarchical clustering.
```
$ nohup Rscript 01_remove_non_redundant_gene_sets.R &
```

# STEP 2: Actually filter based on a dissimilarity value and write a filtered set using this interactive script. We chose a dissimilarity threshold of 0.7. 
```
- 02_examine_gene_set_hierarchical_clusters.ipynb
```
This will produce a GMT file called `filtered_GO_Hs_symbols.gmt` with the filtered gene sets. The 
