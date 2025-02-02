# Run the standard scanpy preprocessing without using scVI integration
# Try harmony vs. no donor-level integration

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


### PART 1: run preprocessing with Harmony integration
start_harmony_time = time.time()

adata.layers['counts'] = adata.X
adata = preprocess_adata(adata, donor_key='donor_id', leiden_resolution=0.5, use_harmony=True)
adata.obsm['X_umap_harmony'] = adata.obsm['X_umap'].copy()

end_harmony_time = time.time()
harmony_elapsed_time = end_harmony_time - start_harmony_time
print(f"Elapsed time for Harmony + other integration steps was {harmony_elapsed_time}")

# perform the leiden clustering; saving the key to "leiden_harmony"
harmony_leiden_key="leiden_harmony"
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution = 0.5, key_added=harmony_leiden_key)

def subsample_per_leiden_cluster(adata, leiden_key, max_cells=1000):

    # create a list to store subsampled indices
    indices_to_keep = []

    # iterate over each Leiden cluster
    for cluster in adata.obs[leiden_key].unique():
        cluster_indices = np.where(adata.obs[leiden_key] == cluster)[0]

        # subsample if cluster size exceeds max_cells
        if len(cluster_indices) > max_cells:
            sampled_indices = np.random.choice(cluster_indices, max_cells, replace=False)
        else:
            sampled_indices = cluster_indices
        indices_to_keep.extend(sampled_indices)

    # Subset the AnnData object
    adata_subsampled = adata[indices_to_keep].copy()
    return adata_subsampled

# subsample to quickly call the marker genes
adata_subsampled = subsample_per_leiden_cluster(adata, max_cells=1000, leiden_key=harmony_leiden_key)
sc.tl.rank_genes_groups(adata_subsampled, groupby=harmony_leiden_key, method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata_subsampled, groupby=harmony_leiden_key, standard_scale="var", n_genes=3)

# transfer the marker genes over to the adata
adata.uns['rank_genes_groups'] = adata_subsampled.uns['rank_genes_groups']
adata.uns['dendrogram_leiden_harmony'] = adata_subsampled.uns['dendrogram_leiden_harmony']




### PART 2: run preprocessing without Harmony integration


print("Performing preprocessing without Harmony integration...", flush=True)
adata.X = adata.layers['counts'] # reset with the raw counts
start_time_no_integration = time.time()
adata = preprocess_adata(adata, donor_key='donor_id', leiden_resolution=0.5, use_harmony=False)
end_time_no_integration = time.time()
elapsed_time_no_integration = end_time_no_integration - start_time_no_integration
print(f"Elapsed time for preprocessing without Harmony is: {elapsed_time_no_integration}", flush=True)
adata.obsm['X_umap_no_harmony'] = adata.obsm['X_umap'].copy()

# save the adata
adata.write("04_harmony_integrated_adata.h5ad")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Entire script completed in {elapsed_time} s!")
