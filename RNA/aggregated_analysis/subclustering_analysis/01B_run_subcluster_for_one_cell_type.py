# perform subclustering for each cell type
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from collections import Counter
import os 
import gc
import time
import scvi

# use argparse with -a <adata_path> -c <cell_type> -o <outdir>
start_time=time.time()

parser = argparse.ArgumentParser(description="Perform subcluster for adata file")
parser.add_argument("-a", "--adata_path", required=True, help="path to adata file")
parser.add_argument("-c", "--cell_type", required=True, help="cell type to subset")
parser.add_argument("-o", "--outdir", required=True, help="output directory")

args = parser.parse_args()

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

# try several different leiden resolutions
resolutions = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]

adata = sc.read_h5ad(args.adata_path)
cell_type = args.cell_type
outdir = args.outdir

# use all genes
num_genes = adata.shape[1]

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=num_genes,
    layer="counts",
    subset=True,
)

sc.tl.pca(adata)

# run scvi
print("Run scVI and train model...", flush=True)

scvi.model.SCVI.setup_anndata(adata, batch_key="tech_plus_study", categorical_covariate_keys = ["donor_id"], layer="counts")

model = scvi.model.SCVI(adata)
model.train(early_stopping=True,  
        max_epochs=500, 
        early_stopping_patience=20, 
        early_stopping_monitor="elbo_validation", 
        train_size=0.9, check_val_every_n_epoch=1, accelerator='cpu')

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent

# store the mean normalized expression in `adata.layers["scvi_normalized"]`
adata.layers["scvi_normalized"] = model.get_normalized_expression(
    adata,
    library_size=10e4,
    transform_batch="Multiome-v1_ENCODE v4 (Snyder)"
)

# use the scvi normalized counts to perform leiden clustering
adata.X = adata.layers['scvi_normalized']
del adata.layers['scvi_normalized']

# call neighbors
SCVI_LATENT_KEY = "X_scVI"
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata)

# perform leiden clustering at each resolution
for resolution in resolutions:
    leiden_key = "leiden" + str(resolution)

    # run leiden clustering on the full adata
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=resolution, key_added = leiden_key)

    # subsample to quickly call the marker genes
    adata_subsampled = subsample_per_leiden_cluster(adata, leiden_key = leiden_key)
    clusters = adata.obs[leiden_key].unique()

    # identify marker genes if there is at least 2 clusters (otherwise, this don't run this snippet)
    if len(clusters) > 1:

        #adata_subsampled = subsample_per_leiden_cluster(adata, leiden_key = leiden_key)

        sc.tl.rank_genes_groups(adata_subsampled, groupby=leiden_key, method="wilcoxon")
        sc.pl.rank_genes_groups_dotplot(adata_subsampled, groupby=leiden_key, standard_scale="var", n_genes=3)

        # transfer the marker genes over to the adata
        adata.uns['rank_genes_groups' + leiden_key] = adata_subsampled.uns['rank_genes_groups']
        #adata.uns['dendrogram_leiden_scVI'] = adata_subsampled.uns['dendrogram_' + leiden_key]
    else:
        print(f"Only 1 cluster found for {leiden_key}. Skipping rank_genes_groups.")

# save the cell-type subsetted adata
adata.write(outdir + "/" + cell_type + ".h5ad")

gc.collect()

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time for this script for is {elapsed_time}", flush=True)
