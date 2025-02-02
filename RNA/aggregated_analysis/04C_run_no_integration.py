import pandas as pd
import numpy as np
import scanpy as sc 
import scanpy.external as sce
import time

def preprocess_adata(adata, donor_key, leiden_resolution, use_harmony=True):
    '''
    Performs:
    1. library size normalization and log scaling
    2. identification of top 2K highly variable genes,
    3. Principal component analysis
    4. Harmony integration (if use_harmony=True)
    5. Neighbor neighbors computation in Harmony integration embedding or the non-Harmony embedding
    6. Leiden clustering (using the KNN graph)

    Parameters:
    adata (AnnData): adata object before preprocessing
    donor_key: the column in adata.obs that corresponds to the donor information (should be 'donor_id')
    leiden_resolution: resolution for leiden clustering, higher means more clusters will be detected

    Returns:
    adata: Postprocessed adata
    '''
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    
    # run without batch key
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    #sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=donor_key)
    sc.tl.pca(adata)

    # use harmony if True
    if use_harmony:
        sce.pp.harmony_integrate(adata, donor_key)
        sc.pp.neighbors(adata, use_rep = "X_pca_harmony")
    else:
        sc.pp.neighbors(adata, use_rep = "X_pca")

    sc.tl.umap(adata)
    sc.pl.umap(adata, color=donor_key, size=2)

    # for this script, skip the leiden clustering part
    #sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution = 0.5)

    return(adata)

start_time = time.time()

print("Loading in the adata...", flush=True)
adata = sc.read_h5ad("03_combined_all_snRNA.h5ad")

### Run preprocessing without any integration
print("Performing preprocessing any integration...", flush=True)
adata.X = adata.layers['counts'] # reset with the raw counts

start_time_no_integration = time.time()
adata = preprocess_adata(adata, donor_key='donor_id', leiden_resolution=0.5, use_harmony=False)
end_time_no_integration = time.time()

elapsed_time_no_integration = end_time_no_integration - start_time_no_integration

print(f"Elapsed time for preprocessing without Harmony is: {elapsed_time_no_integration}", flush=True)
adata.obsm['X_umap_no_harmony'] = adata.obsm['X_umap'].copy()

# save the adata
adata.write("04C_unintegrated_adata.h5ad")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Entire script completed in {elapsed_time} s!")
