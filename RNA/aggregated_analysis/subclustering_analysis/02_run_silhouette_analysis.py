# iteratively run the silhouette analysis

from sklearn.metrics import silhouette_score
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
import os 

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

# iterate through all files in post_subclustering/
post_clust_dir = "post_subclustering/"

# save dfs in silhouette results
sil_dir = "silhouette_res"
os.makedirs(sil_dir, exist_ok=True)

# list all files in folder
files = [f for f in os.listdir(post_clust_dir) if os.path.isfile(os.path.join(post_clust_dir, f))]

# subsampling threshold for # cells (see below)
threshold_for_subsampling=100000

#leiden_clusters =  ["leiden0.1", "leiden0.25", "leiden0.5", "leiden0.75" , "leiden1.0", "leiden1.25", "leiden1.5", "leiden2.0"]
leiden_clusters = ["leiden0.25", "leiden0.5", "leiden0.75" , "leiden1.0", "leiden1.25", "leiden1.5", "leiden2.0"]

for file_name in files: 

    # load in the adata
    adata = sc.read_h5ad(post_clust_dir + file_name)

    # if it contains more than 50K cells, later subsample per leiden cluster to save time
    num_cells = adata.shape[0]
    
    cell_type = os.path.splitext(os.path.basename(file_name))[0]

    print(cell_type)

    sil_score_list = list()

    for leiden in leiden_clusters:
        clusters = adata.obs[leiden].unique()
        
        if len(clusters) <= 1:
            print(f"Only {len(clusters)} cluster(s) in {leiden}, skipping silhouette computation.")
            sil_score_list.append("NA")
            continue
        
        if num_cells < threshold_for_subsampling:
            score = compute_silhouette(adata, cluster_key=leiden, use_rep="X_pca")
        else:
            subsampled_adata = subsample_per_leiden_cluster(adata, leiden_key=leiden, max_cells=3000)
            score = compute_silhouette(subsampled_adata, cluster_key=leiden, use_rep="X_pca")

        sil_score_list.append(score)  # <-- append for both cases

    sil_score_df = pd.DataFrame({
        'leiden': leiden_clusters,
        'sil_score': sil_score_list
    })
    
    sil_score_df.to_csv(sil_dir + "/" + cell_type + ".csv")
