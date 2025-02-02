### Similar to snRNA-seq, we will examine each individual cell type and perform annotation refinement
### We will perform MAGIC imputation to assist with identifying true nuclei for each cell type

import scanpy as sc
import pandas as pd
import numpy as np
import snapatac2 as snap
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
import os
import time
import gc

start_time = time.time()

adata = sc.read_h5ad("05_cleaned_adata.h5ad")

unique_cell_types = adata.obs.cell_type.unique()

# create output directory for individual cell type gene imputed matrices
output_dir = "pre_indiv_cell_type_adata/"
os.makedirs(output_dir, exist_ok=True)

def perform_MAGIC_imputation(peak_adata, leiden_resolution):
    '''Create gene matrix and perform MAGIC imputation'''
    # create gene matrix from peak matrix
    gene_matrix = snap.pp.make_gene_matrix(peak_adata, snap.genome.hg38)
    sc.pp.filter_genes(gene_matrix, min_counts = 0.5)
    sc.pp.normalize_total(gene_matrix)
    sc.pp.log1p(gene_matrix)

    # perform MAGIC imputation
    sc.external.pp.magic(gene_matrix, solver="approximate")
    sc.pp.neighbors(gene_matrix)
    
    # call leiden clusters
    sc.tl.leiden(gene_matrix, key_added='redo_leiden', resolution = leiden_resolution)
    # transfer UMAP from original peak_adata over
    gene_matrix.obsm["X_umap"] = peak_adata.obsm["X_umap"]
    return(gene_matrix)

def add_umap_coords(adata):
    '''Add UMAP1 and UMAP2 to the obs DataFrame for easier manipulation, as scanpy does not show coordinates'''
    umap_coords = adata.obsm["X_umap"]
    adata.obs["UMAP1"] = umap_coords[:, 0]
    adata.obs["UMAP2"] = umap_coords[:, 1]
    return(adata)

for cell_type_val in unique_cell_types:
    print(cell_type_val, flush=True)
    cell_type_adata = adata[adata.obs.cell_type == cell_type_val].copy()
    cell_type_adata = add_umap_coords(cell_type_adata)

    gene_matrix = perform_MAGIC_imputation(cell_type_adata, leiden_resolution = 0.1)

    output_path = output_dir + cell_type_val + "_gene_imputed.h5ad"
    gene_matrix.write(output_path)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script complete! Start time in seconds is {elapsed_time}")
