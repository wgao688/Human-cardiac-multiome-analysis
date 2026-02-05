### trajectory analysis functions

import scipy
import scanpy as sc
from scipy.cluster.hierarchy import fcluster
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import Counter
from scipy import sparse
from scipy.stats import spearmanr
import statsmodels.api as sm
import gseapy as gp
from gseapy import enrichr
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams["pdf.use14corefonts"] = False
from matplotlib.patches import Rectangle
from matplotlib import colors

def downsample_per_group(adata, group_key="age_disease", max_cells=10000,random_state=0,copy=True):

    '''Downsample adata per group'''

    rng = np.random.default_rng(random_state)
    keep_indices = []

    for group, idx in adata.obs.groupby(group_key).indices.items():
        idx = np.array(idx)

        if len(idx) > max_cells:
            sampled = rng.choice(idx, size=max_cells, replace=False)
        else:
            sampled = idx

        keep_indices.append(sampled)

    keep_indices = np.concatenate(keep_indices)

    adata_ds = adata[keep_indices].copy() if copy else adata[keep_indices]

    print("Downsampling summary:")
    print(adata_ds.obs[group_key].value_counts())

    return adata_ds

def compute_gene_signatures_along_traj(adata, trajectory_score_name="traj_score"):

    '''Compute the expression of all genes along the trajectory'''
    
    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)

    traj = adata.obs[trajectory_score_name].values

    rho = np.zeros(X.shape[1])
    pval = np.zeros(X.shape[1])

    # iterate through all genes
    for i in range(X.shape[1]):
        r, p = spearmanr(X[:, i], traj)
        rho[i] = r
        pval[i] = p

    # return df
    traj_df = pd.DataFrame({"gene": adata.var_names,
                            "spearman_r": rho,
                            "pval": pval}).sort_values("spearman_r", ascending=False)
    return traj_df

def perform_trajectory_analysis(adata, other_group, max_cells_per_group=3000, metadata_key = "age_disease",
                                non_diseased_group="postnatal:ND"):
    '''Calculate the density gradient from non-diseased to another state, return the adata object with the gradient.'''

    # subsample the filtered adata to only the select groups
    select_groups = [non_diseased_group, other_group]
    filt_adata = adata[adata.obs[metadata_key].isin(select_groups)].copy()
    # further subsample
    filt_adata_ds = downsample_per_group(filt_adata, group_key=metadata_key, max_cells=max_cells_per_group, random_state=42)

    # recompute the neighbor graph and UMAP using the X_scVI embedding
    sc.pp.neighbors(filt_adata_ds, use_rep = "X_scVI")
    sc.tl.umap(filt_adata_ds)

    # compute the densities for the two groups
    group_adata_list = list()
    for group in select_groups:
        group_only_adata = filt_adata_ds[filt_adata_ds.obs[metadata_key] == group].copy()
        sc.tl.embedding_density(group_only_adata, basis="umap", 
                                groupby=metadata_key, key_added = f"umap_density_{group}")
        group_adata_list.append(group_only_adata)

    combined_adata = sc.concat(group_adata_list, join = "outer")
    density_cols = [f"umap_density_{non_diseased_group}", f"umap_density_{other_group}"]
    combined_adata.obs[density_cols] = combined_adata.obs[density_cols].fillna(0)

    # obtain the densities of each of group
    dens_ND = combined_adata.obs[f"umap_density_{non_diseased_group}"].values.copy()
    dens_other = combined_adata.obs[f"umap_density_{other_group}"].values.copy()

    # calculate the center (max density)
    X = combined_adata.obsm["X_umap"]
    nd_center  = X[dens_ND.argmax()]
    other_center = X[dens_other.argmax()]

    # normalize the densities 
    dn = (dens_ND - dens_ND.min()) / (dens_ND.max() - dens_ND.min())
    dd = (dens_other - dens_other.min()) / (dens_other.max() - dens_other.min())
    
    # disease progression score
    traj_score = dd - dn
    combined_adata.obs["traj_score"] = traj_score

    return combined_adata

def visualize_gene_trajectory(adata, gene, other_state, trajectory_key = "traj_score"):
    '''Visualize gene trajectory along density gradient'''
    
    y = adata[:, gene].X.toarray().ravel() if scipy.sparse.issparse(adata.X) else adata[:, gene].X.ravel()
    traj = adata.obs[trajectory_key]
    lowess = sm.nonparametric.lowess(y, traj, frac=0.2)
    plt.scatter(traj, y, s=4, alpha=0.3)
    plt.plot(lowess[:,0], lowess[:,1], color="red")
    plt.xlabel(f"ND → {other_state} density trajectory")
    plt.ylabel(gene)
    plt.show()

def hierarchical_clustering_traj(adata, trajectory_df,
    other_state, trajectory_key="traj_score", n_quantiles=20, n_labels=30,
    max_plot_genes=3000, figsize_tuple=(15, 20)):

    # use all genes for trajectory analysis
    all_genes = trajectory_df["gene"].values

    # sort cells by trajectory
    traj = adata.obs[trajectory_key].values
    traj_order = np.argsort(traj)
    traj_sorted = traj[traj_order]

    # build expression matrix (ALL genes)
    expr_list = []
    for gene in all_genes:
        y = adata[:, gene].X
        if scipy.sparse.issparse(y):
            y = y.toarray().ravel()
        else:
            y = np.ravel(y)
        expr_list.append(y[traj_order])

    expr_matrix = np.array(expr_list)

    # quantile smoothing
    quantiles = np.linspace(0, 100, n_quantiles + 1)
    expr_quantiled = []

    for i in range(expr_matrix.shape[0]):
        q_gene = []
        for j in range(n_quantiles):
            lo = np.percentile(traj_sorted, quantiles[j])
            hi = np.percentile(traj_sorted, quantiles[j + 1])
            mask = (traj_sorted >= lo) & (traj_sorted <= hi)

            if np.any(mask):
                q_gene.append(np.mean(expr_matrix[i, mask]))
            else:
                q_gene.append(np.nan)

        expr_quantiled.append(q_gene)

    expr_quantiled = np.array(expr_quantiled)

    # get the z-score
    print("Performing Z-score normalization...")
    scaler = StandardScaler()
    expr_zscore = scaler.fit_transform(expr_quantiled.T).T

    # perform hierarchical clustering
    print("Performing hierarchical clustering...")
    linkage_matrix = linkage(expr_zscore, method="ward")
    gene_order = leaves_list(linkage_matrix)

    expr_zscore_clustered = expr_zscore[gene_order, :]
    all_genes_clustered = all_genes[gene_order]

    # distance-based clusters (ALL genes)
    distance_threshold = 0.25 * np.max(linkage_matrix[:, 2])
    cluster_labels = fcluster(linkage_matrix, t=distance_threshold, criterion="distance")
    cluster_labels_ordered = cluster_labels[gene_order]

    print(f"Distance threshold {distance_threshold:.2f} → {len(np.unique(cluster_labels))} clusters")

    # results for ALL genes
    df_results = pd.DataFrame({
        "gene": all_genes,
        "cluster": cluster_labels,
        "spearman_r": trajectory_df["spearman_r"].values,
        "pval": trajectory_df["pval"].values
    })

    # subset genes to plot
    n_genes_total = len(all_genes_clustered)

    if n_genes_total > max_plot_genes:
        plot_idx = np.linspace(0, n_genes_total - 1, max_plot_genes).astype(int)
    else:
        plot_idx = np.arange(n_genes_total)

    expr_plot = expr_zscore_clustered[plot_idx, :]
    genes_plot = all_genes_clustered[plot_idx]
    cluster_plot = cluster_labels_ordered[plot_idx]

    # y-axis labels
    spaced_indices = np.linspace(0, len(genes_plot) - 1, n_labels).astype(int)
    spaced_labels = [genes_plot[i] for i in spaced_indices]

    # produce plot
    fig = plt.figure(figsize=figsize_tuple)

    # dendrogram (ALL genes structure)
    ax_dendro = fig.add_axes([0.02, 0.1, 0.08, 0.8])
    dendrogram(linkage_matrix, orientation="left", no_labels=True,
               ax=ax_dendro, color_threshold=distance_threshold)
    ax_dendro.set_yticks([])
    ax_dendro.set_ylabel("Hierarchical clustering", fontsize=10)

    # cluster bar (subset)
    ax_cluster = fig.add_axes([0.105, 0.1, 0.015, 0.8])
    ax_cluster.set_xlim(0, 1)
    ax_cluster.set_ylim(0, len(cluster_plot))

    cmap = plt.get_cmap("tab20")
    norm = colors.Normalize(vmin=0, vmax=np.max(cluster_labels))

    for i, val in enumerate(cluster_plot):
        ax_cluster.add_patch(Rectangle(
            (0, i), 1, 1,
            facecolor=cmap(norm(val)),
            edgecolor="none", rasterized=True
        ))

    ax_cluster.axis("off")

    # heatmap (subset)
    ax_heatmap = fig.add_axes([0.13, 0.1, 0.75, 0.8])
    sns.heatmap(expr_plot,
                cmap="RdBu_r",
                center=0,
                yticklabels=False,
                xticklabels=False,
                cbar_kws={"label": "Z-score"},
                ax=ax_heatmap, rasterized=True)

    ax_heatmap.set_yticks(spaced_indices + 0.5)
    ax_heatmap.set_yticklabels(spaced_labels, rotation=0, fontsize=9)
    ax_heatmap.yaxis.tick_right()
    ax_heatmap.tick_params(axis="y", length=5, which="major",
                           right=True, labelright=True)

    ax_heatmap.set_xlabel(f"ND to {other_state} trajectory", fontsize=12)
    ax_heatmap.set_title("Gene expression along trajectory", fontsize=14)

    return fig, df_results

def perform_ora_analysis_per_cluster(gene_clusters_df, other_group, pval_threshold = 0.05, gene_set = 'MSigDB_Hallmark_2020'):
    print("Performing Over-Representation Analysis for each cluster...")

    ora_results = {}
    n_clusters = len(gene_clusters_df['cluster'].unique())
    
    for i in range(1, n_clusters + 1):
        cluster_genes = gene_clusters_df[gene_clusters_df['cluster'] == i]['gene'].tolist()
        print(f"\nCluster {i}: {len(cluster_genes)} genes")
        # perform enrichment analysis
        try:
            enr = gp.enrichr(gene_list=cluster_genes, gene_sets=gene_set, organism='Human',
                outdir=None, cutoff=0.05)
            ora_results[f'Cluster_{i}'] = enr.results
        
        except Exception as e:
            print(f"  Error in enrichment analysis: {e}")

    ora_res_df = pd.concat(ora_results)

    # filter for significant hits
    ora_res_df_sig = ora_res_df[ora_res_df['Adjusted P-value'] < pval_threshold].copy()
    ora_res_df_sig = ora_res_df_sig.reset_index().rename(columns = {'level_0': 'cluster'}).drop(columns = ['level_1'])
    ora_res_df_sig['contrast'] = other_group

    return ora_res_df_sig

def run_trajectory_pipeline(adata, other_group, max_cells_per_group=3000, non_diseased_group = "postnatal:ND", num_quantiles=20):

    '''Run the full trajectory analysis from non-diseased to another group'''

    # get the trajectory values
    print("Calculating the densities...", flush=True)
    comparison_adata = perform_trajectory_analysis(adata, other_group = other_group,
                                               max_cells_per_group=max_cells_per_group, non_diseased_group=non_diseased_group)

    print("Track gene expression along gradient...", flush=True)
    # get the gene expression across the trajectory
    gene_traj_df = compute_gene_signatures_along_traj(comparison_adata)

    # perform hierarchical cluster; produce heatmap
    print("Performing hierarchical clustering...", flush=True)
    fig, df_results = hierarchical_clustering_traj(adata = comparison_adata,
                                                   other_state = other_group, trajectory_df = gene_traj_df, n_quantiles = num_quantiles)
    # perform ora for each cluster
    print("Performing ORA...", flush=True)
    ora_res_df = perform_ora_analysis_per_cluster(df_results, other_group = other_group)

    return df_results, ora_res_df, fig
