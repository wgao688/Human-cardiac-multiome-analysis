# iteratively run the silhouette analysis

from sklearn.metrics import silhouette_score
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
import os 
import time
import argparse

parser = argparse.ArgumentParser(description="Load an AnnData object for subclustering analysis")
parser.add_argument("-c", "--cell_type", type=str, required=True, help="the name of the cell type")

args = parser.parse_args()

start_time = time.time()

def subsample_per_leiden_cluster(adata, leiden_key, max_cells=1000):

    '''
    Subsample for cell types with large number of cells
    '''

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

def compute_silhouette(adata, cluster_key, use_rep="X_pca"):
    '''
    Compute the mean silhouette score
    '''

    # extract embedding and cluster labels
    X = adata.obsm[use_rep]
    labels = adata.obs[cluster_key].values

    score = silhouette_score(X, labels)
    return score

# get the full path 
sil_dir = "silhouette_res_scVI"
os.makedirs(sil_dir, exist_ok=True)

# subsampling threshold for # cells (see below)
threshold_for_subsampling=50000

#leiden_clusters =  ["leiden0.1", "leiden0.25", "leiden0.5", "leiden0.75" , "leiden1.0", "leiden1.25", "leiden1.5", "leiden2.0"]
#leiden_clusters = ["leiden0.25", "leiden0.5", "leiden0.75" , "leiden1.0", "leiden1.25", "leiden1.5", "leiden2.0"]
leiden_clusters = ["leiden0.25", "leiden0.5", "leiden0.75" , "leiden1.0", "leiden1.25", "leiden1.5", "leiden1.75"]

# load in the adata
cell_type = args.cell_type
adata_path = "continuous_cov_corrected_no_scvi/" + args.cell_type + ".h5ad"
adata = sc.read_h5ad(adata_path)

# if it contains more than 50K cells, later subsample per leiden cluster to save time
num_cells = adata.shape[0]
    
#cell_type = os.path.splitext(os.path.basename(file_name))[0]

print(cell_type)

sil_score_list = list()

for leiden in leiden_clusters:
    clusters = adata.obs[leiden].unique()
        
    if len(clusters) <= 1:
        print(f"Only {len(clusters)} cluster(s) in {leiden}, skipping silhouette computation.")
        sil_score_list.append("NA")
        continue
        
    if num_cells < threshold_for_subsampling:
        #score = compute_silhouette(adata, cluster_key=leiden, use_rep="X_pca")
        score = compute_silhouette(adata, cluster_key=leiden, use_rep="X_scVI")
    else:
        subsampled_adata = subsample_per_leiden_cluster(adata, leiden_key=leiden, max_cells=1000)
        #score = compute_silhouette(subsampled_adata, cluster_key=leiden, use_rep="X_pca")
        score = compute_silhouette(subsampled_adata, cluster_key=leiden, use_rep="X_scVI")

    sil_score_list.append(score)  # append for both cases

sil_score_df = pd.DataFrame({
    'leiden': leiden_clusters,
    'sil_score': sil_score_list
})
    
sil_score_df.to_csv(sil_dir + "/" + cell_type + ".csv")

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Elapsed time for script is {elapsed_time} s", flush=True)
