# Performs gene matrix imputation on the harmonized adata

import snapatac2 as snap
import scanpy as sc
import numpy as np
import tempfile
import os
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time

start_time = time.time()

print("Reading in the harmonized adata...", flush=True)
adata=sc.read_h5ad("harmonized_adata.h5ad")

### PERFORM GENE MATRIX IMPUTATION 
print("starting making gene matrix...", flush=True)
# create gene matrix
gene_matrix = snap.pp.make_gene_matrix(adata, snap.genome.hg38)
print("done making gene matrix")
# save the raw counts
gene_matrix.layers['raw_counts'] = gene_matrix.X

mid_time = time.time()

# copy over the UMAP from the spectral embedding
gene_matrix.obsm["X_umap"] = adata.obsm["X_umap"]

# perform differential gene expression analysis
sc.tl.rank_genes_groups(gene_matrix, groupby="leiden", method="wilcoxon")

# save the imputed adata
gene_matrix.write("gene_matrix_non_MAGIC_imputed.h5ad")

elapsed_time_1 = mid_time - start_time
print(f"Elapsed time for gene imputation without MAGIC is {elapsed_time_1}", flush=True)
print("Continue with MAGIC imputation...", flush=True)

#### MAGIC IMPUTATION ADDITIONAL STEPS
# perform MAGIC imputation, which requires log transformed data
sc.pp.filter_genes(gene_matrix, min_cells = 1000) # remove lowly expressed genes
sc.pp.normalize_total(gene_matrix)
sc.pp.log1p(gene_matrix)

# this may take quite a while
# For extremely large datasets, using n_pca < 20 allows neighborhoods to be calculated in roughly log(n_samples) time.
sc.external.pp.magic(gene_matrix, solver="approximate", n_pca = 15) # use n_pca = 15 since this is a large dataset;
# “approximate” uses a faster implementation that performs imputation in the PCA space and then projects back to the gene space.
# perform differential gene expression analysis with MAGIC imputation
sc.tl.rank_genes_groups(gene_matrix, groupby="leiden", method="wilcoxon")

# save the imputed adata
gene_matrix.write("gene_matrix_MAGIC_imputed.h5ad")

end_time = time.time()
elapsed_time_2 = end_time - start_time

print(f"Elapsed time for the whole script is {elapsed_time_2}. Script complete!", flush=True)
