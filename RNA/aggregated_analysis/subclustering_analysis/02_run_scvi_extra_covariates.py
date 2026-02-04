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

def calculate_covariates(adata, layer='counts'):
    '''Calculate covariates that we want to regress out'''

    # calculate log10_UMI
    adata.obs['total_UMI'] = adata.layers[layer].sum(axis = 1)
    adata.obs['log10_UMI'] = np.log10(adata.obs['total_UMI'])

    # calculate Endothelial, CM, and Fib genes
    endo_genes = ["PECAM1", "VWF", "KDR", "FLT1", "ESAM", "EMCN", "RAMP2", "CD34", "ENG", "PLVAP"]
    cm_genes = ["TNNT2", "TTN", "MYH6", "MYH7", "ACTC1", "TNNI3", "MYL2", "MYL7", "RYR2", "PLN"]
    fib_genes = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL5A1", "COL6A1", "COL6A2", "COL6A3", "FBLN1"]
    endo_genes = [g for g in endo_genes if g in adata.var_names]
    cm_genes   = [g for g in cm_genes   if g in adata.var_names]
    fib_genes = [g for g in fib_genes if g in adata.var_names]

    sc.tl.score_genes(adata, endo_genes, score_name="endo_score")
    sc.tl.score_genes(adata, cm_genes, score_name="cm_score")
    sc.tl.score_genes(adata, fib_genes, score_name="fib_score")

    # calculate pct_mt, ribo, hb
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    return(adata)

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

os.makedirs(outdir, exist_ok=True)

# use all genes
num_genes = adata.shape[1]

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=num_genes, layer="counts", subset=True)

# calculate the additional covariates
adata = calculate_covariates(adata)
sc.tl.pca(adata)

# build continuous covariates list
base_covs = ["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb", "log10_UMI"]
lineage_covs = ["endo_score", "fib_score", "cm_score"]

if cell_type == "Endothelial":
    continuous_covariates = ["fib_score", "cm_score"] + base_covs

elif cell_type == "Fibroblast":
    continuous_covariates = ["endo_score", "cm_score"] + base_covs

elif cell_type == "Cardiomyocyte":
    continuous_covariates = ["endo_score", "fib_score"] + base_covs

else:
    # all other cell types: regress out all three lineage scores
    continuous_covariates = lineage_covs + base_covs

# run scvi
print("Run scVI and train model...", flush=True)
scvi.model.SCVI.setup_anndata(adata, batch_key="tech_plus_study", categorical_covariate_keys = ["donor_id"], 
        continuous_covariate_keys = continuous_covariates,
        layer="counts")

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

    # identify differential genes across the leiden clusters
    #sc.tl.rank_genes_groups(subset_adata, groupby=leiden_key, method="wilcoxon", key_added = leiden_key)

# save the cell-type subsetted adata
adata.write(outdir + "/" + cell_type + ".h5ad")

gc.collect()

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time for this script for is {elapsed_time}", flush=True)
