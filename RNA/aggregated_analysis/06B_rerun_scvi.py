import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from collections import Counter
import os 
import gc
import time
import scvi

start_time = time.time()

print("Loading in the adata...", flush=True)
adata = sc.read_h5ad("06A_filtered_adata.h5ad")

# use all genes
num_genes = adata.shape[1]
#sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=num_genes, layer="counts", subset=True)
#sc.tl.pca(adata)

# build continuous covariates list
base_covs = ["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb", "log10_UMI"]
lineage_covs = ["endo_score", "fib_score", "cm_score"]
continuous_covariates = lineage_covs + base_covs

# run scvi
print("Run scVI and train model...", flush=True)
scvi.model.SCVI.setup_anndata(adata, batch_key="tech_plus_study", categorical_covariate_keys = ["donor_id"], 
        continuous_covariate_keys = continuous_covariates,
        layer="counts")

model = scvi.model.SCVI(adata)
model.train(early_stopping=True, early_stopping_monitor="elbo_validation",  check_val_every_n_epoch=1, accelerator='cpu')

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent

# recompute neighbors using scVI latent
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)

# recompute UMAP
sc.tl.umap(adata)

# store normalized in adata.X
adata.X = model.get_normalized_expression(
    adata,
    library_size=10e4,
    transform_batch="Multiome-v1_ENCODE v4 (Snyder)"
)

# save the cell-type subsetted adata
adata.write("06B_combined_scVI.h5ad")

# save the one without the scVI layer
del adata.X
adata.write("06B_combined_scVI_raw_only.h5ad")

gc.collect()

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time for this script for is {elapsed_time}", flush=True)
