import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import gseapy as gp
import os
import pickle


def produce_crosstab_barplot(metadata, cluster_label, metadata_column):
    '''Produce a barplot showing how the cluster_label (e.g. leiden) distributes according to a metadata column (e.g. study)'''

    prop_df = pd.crosstab(metadata[cluster_label], metadata[metadata_column])
    norm_prop_df = prop_df.div(prop_df.sum(axis = 1), axis = 0)

    ax = norm_prop_df.plot(kind = 'bar', stacked=True, figsize=(8, 4), colormap='tab20')
    ax.set_xlabel(cluster_label)
    ax.set_ylabel("Proportion")
    ax.set_title(f"Proportion per {cluster_label} grouped by {metadata_column}")

    # reverse the legend order so that it matches up with the barplot
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title={cluster_label}, bbox_to_anchor=(1.05, 1), loc='upper left')

    # rotate the x labels
    plt.xticks(rotation=90, ha='right')
    return plt

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


def load_cnMF(cell_type):

    parent_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))

    file_path = os.path.join(parent_dir, "cNMF", cell_type + cell_type)
    
    with open(file_path, "rb") as f:
        cnmf_obj = pickle.load(f)

    cnmf_obj.k_selection_plot()

    return cnmf_obj


def silhouette(adata, cell_type, cell_threshold=50000):
    '''Plot the silhouette analysis'''

    # load in the silhouette score results
    sil_score_df = pd.read_csv("../silhouette_res/" + cell_type + ".csv", index_col = 0)

    sil_plot = sns.scatterplot(sil_score_df, x = 'leiden', y = 'sil_score')
    
    return(sil_plot)

def specify_leiden_subclustering(adata, leiden, cell_threshold=50000):

    '''Identify the best leiden resolution. Subsample for marker gene analysis if above cell_threshold for faster computation'''

    top_resolution = leiden

    # produce barplot with sil score and study
    produce_crosstab_barplot(adata.obs, cluster_label=top_resolution, metadata_column = "study")
    produce_crosstab_barplot(adata.obs, cluster_label=top_resolution, metadata_column = "age_status")
    produce_crosstab_barplot(adata.obs, cluster_label=top_resolution, metadata_column = "disease_binary")

    sc.pl.umap(adata, color = ["study", top_resolution], legend_loc = "on data")

    num_cells = adata.shape[0]

    # get marker genes for the best resolution
    if num_cells > cell_threshold:
        subsampled_adata = subsample_per_leiden_cluster(adata, leiden_key = top_resolution, max_cells = 5000)
        sc.tl.rank_genes_groups(subsampled_adata, groupby=top_resolution, method="wilcoxon")
        sc.tl.dendrogram(subsampled_adata, groupby=top_resolution)
        sc.pl.rank_genes_groups_dotplot(subsampled_adata, groupby=top_resolution, standard_scale="var", n_genes=5)
        adata.uns['rank_genes_groups'] = subsampled_adata.uns['rank_genes_groups']
    else:
        sc.tl.rank_genes_groups(adata, groupby=top_resolution, method="wilcoxon")
        sc.tl.dendrogram(adata, groupby=top_resolution)
        sc.pl.rank_genes_groups_dotplot(adata, groupby=top_resolution, standard_scale="var", n_genes=5)

    print(Counter(adata.obs[top_resolution]))

    return [top_resolution, adata]

def produce_cNMF_leiden_overlap(adata, NMF_obj, k_opt, leiden_key):
    NMF_obj.consensus(k=k_opt, density_threshold=0.01)

    usage, spectra_scores, spectra_tpm, top_genes = NMF_obj.load_results(
        K=k_opt, density_threshold=0.01
    )

    subsampled_adata = adata[usage.index].copy()

    subsampled_adata.obsm["X_cNMF"] = usage.loc[subsampled_adata.obs_names].values
    subsampled_adata.uns["cNMF_programs"] = spectra_scores 
    subsampled_adata.uns["cNMF_top_genes"] = top_genes

    # produce heatmap of mean program values 
    programs_df = pd.DataFrame(
        subsampled_adata.obsm["X_cNMF"], 
        index=subsampled_adata.obs_names, 
        columns=[str(i+1) for i in range(subsampled_adata.obsm["X_cNMF"].shape[1])]
    )

    # add Leiden cluster info
    programs_df["leiden"] = subsampled_adata.obs[leiden_key]

    # compute cluster means
    cluster_means = programs_df.groupby("leiden").mean()

    # heatmap
    plt.figure(figsize=(8,5))
    sns.heatmap(cluster_means, cmap="viridis")
    plt.title("Average cNMF program activity per Leiden cluster")
    plt.xlabel("Program")
    plt.ylabel("Leiden cluster")

    return plt, top_genes, programs_df
