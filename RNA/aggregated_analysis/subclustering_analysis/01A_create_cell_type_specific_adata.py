# generate individual cell type adata

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from collections import Counter
import os 
import gc
import time
import scvi

pre_clustering_dir = "pre_subclustering/"
os.makedirs(pre_clustering_dir, exist_ok=True)

post_clustering_dir = "post_subclustering/"
os.makedirs(post_clustering_dir, exist_ok=True)

adata_path = "../05E_all_snRNA_adata_without_scvi.h5ad"
print("Loading in the annotated adata...", flush=True)
adata = sc.read_h5ad(adata_path)

cell_types_to_subcluster = adata.obs.final_cell_type.unique()

cell_type_key="v2_scvi_cell_type"

# save each adata before subclustering
for cell_type in cell_types_to_subcluster:

    print(cell_type, flush=True)

    subset_adata = adata[adata.obs[cell_type_key] == cell_type].copy()

    subset_adata_path = pre_clustering_dir + "/" + cell_type + ".h5ad"
    subset_adata.write(subset_adata_path)
