import scanpy as sc
import scvi
import time
import numpy as np
import pandas as pd

# Perform scVI integration

start_time = time.time()

# import the adata
print("Loading in the adata...", flush=True)
adata = sc.read_h5ad("03_combined_all_snRNA.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#num_genes = 2000
# use all the genes
num_genes = adata.shape[1]

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=num_genes,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
    )

scvi_start_time = time.time()

print("Specifying and training the model...", flush=True)
scvi.model.SCVI.setup_anndata(
    adata, batch_key="tech_plus_study",
    categorical_covariate_keys = ["donor_id"],
    layer="counts",
)

model = scvi.model.SCVI(adata)
model.train()
model.save("scVI_model/", overwrite = True)

scvi_end_time = time.time()
scvi_elapsed_time = scvi_end_time - scvi_start_time
print(f"Elapsed time for scvi integration was {scvi_elapsed_time}", flush=True)

# get the latent embedding and add this to the adata object
latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent

# Store the mean normalized expression in `adata.layers["scvi_normalized"]`
adata.layers["scvi_normalized"] = model.get_normalized_expression(
        adata,
        library_size=10e4,
        transform_batch="Multiome-v1_ENCODE v4 (Snyder)"
    )

batch_col = "donor_id"
corrected_data = model.get_normalized_expression(transform_batch = sorted(adata.obs[batch_col].unique()),
                                                 library_size = 1e4)

# use scVI latent space for computing KNN graph and UMAP generation
print("Performing KNN and UMAP...",flush=True)
SCVI_LATENT_KEY = "X_scVI"
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=0.3)

print("Performing leiden clustering...", flush=True)
SCVI_CLUSTERS_KEY = "leiden_scVI"
sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=0.5, flavor="igraph", n_iterations=2)
sc.pl.umap(adata, color=[SCVI_CLUSTERS_KEY], frameon=False)

# Run the following in a subsampled manner per leiden cluster in a later script to reduce amount of time 
# Obtain cluster-specific differentially expressed genes
print("Identifying marker genes for each cluster...", flush=True)

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
adata_subsampled = subsample_per_leiden_cluster(adata, max_cells=1000, leiden_key=SCVI_CLUSTERS_KEY)
sc.tl.rank_genes_groups(adata_subsampled, groupby=SCVI_CLUSTERS_KEY, method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata_subsampled, groupby=SCVI_CLUSTERS_KEY, standard_scale="var", n_genes=3)

# transfer the marker genes over to the adata
adata.uns['rank_genes_groups'] = adata_subsampled.uns['rank_genes_groups']
adata.uns['dendrogram_leiden_scVI'] = adata_subsampled.uns['dendrogram_leiden_scVI']

adata.write("05_scvi_adata.h5ad")

end_time = time.time()
elapsed_time = end_time - start_time

print("Elapsed time for script in seconds is: " + str(elapsed_time), flush=True)
print("Script complete!", flush=True)
