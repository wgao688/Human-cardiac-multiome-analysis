# Attempt scvi integration of snRNA-seq and scRNA-seq datasets

import scanpy as sc
import scvi
from collections import Counter
import time
import numpy as np
import pandas as pd
sc._settings.settings._vector_friendly=True
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

# Perform scVI integration
start_time = time.time()


adata = sc.read_h5ad("03B_sc_sn_adata.h5ad")

# preprocess and setup scVI object
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

num_genes = adata.shape[1]

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=num_genes,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
    )

# make sure to include cell_or_nuclei as a categorical covariate
scvi.model.SCVI.setup_anndata(
    adata, batch_key="tech_plus_study",
    categorical_covariate_keys = ["donor_id", "cell_or_nuclei"],
    layer="counts",
)

model = scvi.model.SCVI(adata)
model.train()
model.save("scVI_cell_v_nuclei/", overwrite = True)

# get the latent embedding and add this to the adata object
latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent

# Store the mean normalized expression in `adata.layers["scvi_normalized"]`
adata.layers["scvi_normalized"] = model.get_normalized_expression(
        adata,
        library_size=10e4,
        transform_batch="5prime-v1_Koenig 2022"
    )

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

adata.write("04A_scvi_sc_sn_combined.h5ad")

end_time = time.time()
elapsed_time = end_time - start_time

plots_dir="scvi_plots/"

# produce UMAP plots
with plt.rc_context():
    sc.pl.umap(adata, color=["cell_or_nuclei"], frameon=False, show = False)
    plt.savefig(plots_dir + "UMAP_cell_vs_nuclei.pdf", bbox_inches="tight")
    plt.show()

with plt.rc_context():
    sc.pl.umap(adata, color=["study"], frameon=False, show = False)
    plt.savefig(plots_dir + "UMAP_study.pdf", bbox_inches="tight")
    plt.show()

with plt.rc_context():
    sc.pl.umap(adata, color=["consistent_cell_type"], frameon=False, show = False, legend_loc = "on data")
    plt.savefig(plots_dir + "UMAP_cell_type.pdf", bbox_inches="tight")
    plt.show()

with plt.rc_context():
    sc.pl.umap(adata, color=["technology"], frameon=False, show = False)
    plt.savefig(plots_dir + "UMAP_technology.pdf", bbox_inches="tight")
    plt.show()

print("Elapsed time for script in seconds is: " + str(elapsed_time), flush=True)
print("Script complete!", flush=True)
