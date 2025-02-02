# Send the subclustering for each cell type

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from collections import Counter
import os 
import gc
import time

start_time=time.time()

cell_type_adata_dir = "scvi_pre_cleanup_subset_cell_type_adata/"
os.makedirs(cell_type_adata_dir, exist_ok=True)

print("Loading in the annotated adata...", flush=True)
adata = sc.read_h5ad("05B_filtered_scvi.h5ad")

# run subclustering for each of the cell types
cell_types_to_subcluster = adata.obs.scvi_cell_type.unique()

# iterate through each of the cell types are perform cell type clustering
for cell_type in cell_types_to_subcluster:

    print(cell_type, flush=True)
    
    # subset the adata to that cell type
    subset_adata = adata[adata.obs['scvi_cell_type'] == cell_type, :]

    # use the scvi normalized counts to perform leiden clustering
    subset_adata.layers['scvi_normalized'] = subset_adata.X

    # leiden clustering - use resolution = 0.5
    sc.tl.leiden(subset_adata, flavor="igraph", n_iterations=2, resolution=0.5, key_added = "redo_leiden_0.5")

    # identify differential genes across the leiden clusters
    sc.tl.rank_genes_groups(subset_adata, groupby="redo_leiden_0.5", method="wilcoxon")

    # save the subsetted adata
    subset_adata.write(cell_type_adata_dir + cell_type + ".h5ad")

    gc.collect()

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time for this script is {elapsed_time}", flush=True)
