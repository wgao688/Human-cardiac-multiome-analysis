import scanpy as sc
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import gseapy as gp
import os
import pickle
from scipy.stats import chisquare
from statsmodels.stats.multitest import multipletests

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

    # subset the adata
    adata_subsampled = adata[indices_to_keep].copy()
    return adata_subsampled

def load_cnMF(cell_type):

    parent_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))

    file_path = os.path.join(parent_dir, "cNMF", cell_type + cell_type)
    
    with open(file_path, "rb") as f:
        cnmf_obj = pickle.load(f)

    cnmf_obj.k_selection_plot()

    return cnmf_obj


def silhouette(adata, cell_type, silhouette_res_dir, cell_threshold=50000):
    '''Plot the silhouette score analysis results'''

    # load in the silhouette score results
    sil_score_df = pd.read_csv(silhouette_res_dir + cell_type + ".csv", index_col = 0)

    plt.figure(figsize=(8,5))
    sns.scatterplot(data=sil_score_df, x='leiden', y='sil_score', s=50, color='blue', alpha=0.7)
    sns.lineplot(data=sil_score_df, x='leiden', y='sil_score', color='red', marker='o')

    plt.xlabel("Leiden Cluster")
    plt.ylabel("Silhouette Score")
    plt.title(f"Silhouette Analysis for {cell_type}")
    plt.xticks(rotation=45)  # rotate x-axis labels
    plt.tight_layout()

    return plt

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


def leiden_resolution_summary(adata, make_plot=True, figsize=(6,4)):
    '''Count number of clusters across all leiden resolutions in adata.obs. Returns df and plot'''

    # find all leiden clusters such as "leiden0.1" except for "leiden_scVI"
    leiden_cols = [c for c in adata.obs.columns if c.startswith("leiden") and "scvi" not in c.lower()]
    
    if len(leiden_cols) == 0:
        raise ValueError("No columns starting with 'leiden' found in adata.obs")

    # count number of unique clusters
    cluster_counts = {col: adata.obs[col].nunique() for col in leiden_cols}

    cluster_df = pd.DataFrame.from_dict( cluster_counts, orient="index", columns=["n_clusters"]).reset_index().rename(columns={"index": "resolution"})

    # sort by cluster number
    cluster_df["res_numeric"] = (
        cluster_df["resolution"]
        .str.replace("leiden", "", regex=False)
        .astype(float)
    )
    cluster_df = cluster_df.sort_values("res_numeric")

    fig = None

    if make_plot:
        fig, ax = plt.subplots(figsize=figsize)

        sns.lineplot(data=cluster_df, x="res_numeric", y="n_clusters", marker="o", ax=ax)

        ax.set_xlabel("Leiden resolution")
        ax.set_ylabel("Number of clusters")
        ax.set_title("Number of clusters across Leiden resolutions")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

    return cluster_df, fig

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

def produce_crosstab_barplot(metadata, cluster_label, metadata_column, figsize_tuple = (8, 4)):
    '''Produce a barplot showing how the cluster_label (e.g. leiden) distributes according to a metadata column (e.g. study)'''

    prop_df = pd.crosstab(metadata[cluster_label], metadata[metadata_column])
    norm_prop_df = prop_df.div(prop_df.sum(axis = 1), axis = 0)

    ax = norm_prop_df.plot(kind = 'bar', stacked=True, figsize=figsize_tuple, colormap='tab20')
    ax.set_xlabel(cluster_label)
    ax.set_ylabel("Proportion")
    ax.set_title(f"Proportion per {cluster_label} grouped by {metadata_column}")

    # reverse the legend order so that it matches up with the barplot
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title={cluster_label}, bbox_to_anchor=(1.05, 1), loc='upper left')

    # rotate the x labels
    plt.xticks(rotation=90, ha='right')
    return plt

def run_enrichr_ora(adata, group, top_n, gene_sets = ["MSigDB_Hallmark_2020"], organism = "Human", outdir = None, print_top = 5):
    '''Perform ORA (Over-Representation Analysis) on the top marker genes of a scanpy cluster using Enrichr.'''
    
    # get top N marker genes
    genes = sc.get.rank_genes_groups_df(adata, group=group).head(top_n)["names"].str.upper().tolist()
    
    # run enrichr ORA
    enr = gp.enrichr(gene_list=genes, gene_sets=gene_sets, organism=organism, outdir=outdir, cutoff=0.05)
    
    # print top results
    #if print_top > 0:
    #    print(enr.results.head(print_top))
    
    return enr

def plot_umap_by_cluster(adata, condition_key, max_per_condition=10000, random_state=42, min_cells=50):
    '''Plot UMAP density maps for each group in condition_key, with downsampling per group.'''
    # downsample per condition
    rng = np.random.default_rng(random_state)
    keep_indices = []

    for group, idx in adata.obs.groupby(condition_key).indices.items():
        idx = np.array(idx)

        # skip very small groups
        if len(idx) < min_cells:
            print(f"Skipping {group}: only {len(idx)} cells")
            continue

        if len(idx) > max_per_condition:
            idx = rng.choice(idx, size=max_per_condition, replace=False)

        # variance check in UMAP space
        umap = adata.obsm["X_umap"][idx]
        if np.var(umap[:, 0]) == 0 or np.var(umap[:, 1]) == 0:
            print(f"Skipping {group}: zero UMAP variance")
            continue

        keep_indices.append(idx)

    keep_indices = np.concatenate(keep_indices)
    adata_ds = adata[keep_indices]
    adata_ds.X = None
    adata_ds.layers.clear()

    print("Downsampling summary:")
    print(adata_ds.obs[condition_key].value_counts())

    # produce density plot
    sc.tl.embedding_density(adata_ds, basis="umap", groupby = condition_key)
    sc.pl.embedding_density(adata_ds, basis="umap", key=f"umap_density_{condition_key}")

def plot_all_cluster_enrichment(adata, cluster_key='leiden', outdir="enrichment_results", **kwargs):

    '''Produce plot for GSEA analysis across all clusters'''

    all_results = []
    
    # sort clusters
    clusters = sorted(adata.obs[cluster_key].unique().tolist(), key=lambda x: int(x))
    
    for cluster in clusters:
        print(f"Processing Cluster {cluster}...", end=" ")
        
        try:
            # run enrichment with error handling for the "No enrich terms" 
            enr_obj = run_enrichr_ora(adata, group=str(cluster), outdir=outdir, **kwargs)
            
            # extract results if they exist
            if enr_obj is not None and hasattr(enr_obj, 'results'):
                res_df = enr_obj.results
                if not res_df.empty:
                    res_df['cluster'] = str(cluster)
                    all_results.append(res_df)
                    print("Success.")
                else:
                    print("Empty results.")
            else:
                print("No results object.")
                
        except ValueError as e:
            # catch "No enrich terms when cutoff = 0.05" error specifically
            print(f"Skipped (No significant terms found).")
            continue
        except Exception as e:
            print(f"Skipped due to unexpected error: {e}")
            continue

    # combine and plot
    if not all_results:
        print("\nRESULT: No significant hits found across ANY clusters. Plotting aborted.")
        return
        
    df_summary = pd.concat(all_results)
    
    # filter on significance, unless this is all empty
    df_summary = df_summary[df_summary['Adjusted P-value'] < 0.05].copy()
    if df_summary.empty:
        print("\nRESULT: No terms passed Adj P-value < 0.05. Plotting aborted.")
        return

    # visualize the top 5 hits per cluster for the summary dotplot
    top_terms = df_summary.groupby('cluster').head(5)['Term'].unique()
    plot_df = df_summary[df_summary['Term'].isin(top_terms)].copy()
    plot_df['nlog10'] = -np.log10(plot_df['Adjusted P-value'].replace(0, 1e-10))

    # create dotplot
    plt.figure(figsize=(10, len(top_terms) * 0.3 + 2))
    sns.scatterplot(data=plot_df, x='cluster', y='Term', size='nlog10',
        hue='Combined Score', palette='viridis', sizes=(40, 300), edgecolor='black')
    plt.title("Summary Enrichment Dotplot", fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="-log10(Adj P)")
    plt.tight_layout()
    
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        plt.savefig(f"{outdir}/combined_dotplot.png", dpi=300)
    plt.show()


def chi_square_per_subcluster(metadata, figsize_tuple, cluster_label, study_col,
                              subcluster_order=None, metadata_order=None,
                              min_cells=30, qval_cutoff=0.05,
                              fc_star_thresh=0.5, plot=True):

    # get the global counts and total cells
    global_counts = metadata[study_col].value_counts().sort_index()
    N = len(metadata)

    test_records = []
    enrich_records = []

    # loop over clusters
    for clust, sub in metadata.groupby(cluster_label):
        if len(sub) < min_cells:
            continue

        obs = sub[study_col].value_counts().reindex(global_counts.index, fill_value=0)
        exp = global_counts * (len(sub) / N)

        # chi-square test for the whole cluster
        stat, p = chisquare(f_obs=obs, f_exp=exp)

        test_records.append({
            "subcluster": clust,
            "chi2": stat,
            "pvalue": p,
            "n_cells": len(sub)
        })

        # loop over metadata categories and skip if below min_cells
        for cat in global_counts.index:

            # if there aren't enough cells (specified by min_cells; then set FC to 0)
            if obs[cat] < min_cells:
                log2fc = 0
                odds_ratio = 1.0

            enrich_records.append({
                "subcluster": clust,
                "metadata_category": cat,
                "observed": obs[cat],
                "expected": exp[cat],
                "odds_ratio": (obs[cat] + 1e-6) / (exp[cat] + 1e-6),
                "log2_fc": np.log2((obs[cat] + 1e-6) / (exp[cat] + 1e-6)),
                "n_in_subcluster": len(sub)
            })

    # compile enrichment results
    test_df = pd.DataFrame(test_records)
    test_df["padj"] = multipletests(test_df["pvalue"], method="fdr_bh")[1]

    enrich_df = pd.DataFrame(enrich_records)
    enrich_df = enrich_df.merge(
        test_df[["subcluster", "padj"]],
        on="subcluster",
        how="left"
    )

    if plot:
        plot_df = enrich_df.query("n_in_subcluster >= @min_cells").copy()
        plot_df["subcluster"] = plot_df["subcluster"].astype(str)
        plot_df["metadata_category"] = plot_df["metadata_category"].astype(str)

        # enforce metadata (column) order if provided
        if metadata_order is not None:
            metadata_order = [str(x) for x in metadata_order]
            keep = [x for x in metadata_order if x in plot_df["metadata_category"].unique()]
            plot_df["metadata_category"] = pd.Categorical(
                plot_df["metadata_category"],
                categories=keep,
                ordered=True
            )
        else:
            keep = sorted(plot_df["metadata_category"].unique())

        # enforce subcluster (row) order
        if subcluster_order is not None:
            subcluster_order = [str(x) for x in subcluster_order]
            leiden_order = [x for x in subcluster_order if x in plot_df["subcluster"].unique()]
        else:
            try:
                leiden_order = sorted(plot_df["subcluster"].astype(int).unique())
                leiden_order = [str(x) for x in leiden_order]
            except ValueError:
                leiden_order = sorted(plot_df["subcluster"].unique())

        plot_df["subcluster"] = pd.Categorical(
            plot_df["subcluster"],
            categories=leiden_order,
            ordered=True
        )

        # pivot wide for heatmap
        mat = plot_df.pivot(index="subcluster", columns="metadata_category", values="log2_fc")
        qmat = plot_df.pivot(index="subcluster", columns="metadata_category", values="padj")

        # enforce column order
        mat = mat.reindex(columns=keep)
        qmat = qmat.reindex(columns=keep)

        mask = qmat >= qval_cutoff

        plt.figure(figsize=figsize_tuple)

        # make the colormap consistent across all samples
        ax = sns.heatmap(
            mat, cmap="coolwarm", center=0, vmin=-2, vmax=2,
            mask=mask, linewidths=0.5, linecolor="lightgrey",
            cbar_kws={"label": "log2(observed / expected)"}
        )

        # add stars for strong effects
        #for i in range(mat.shape[0]):
        #    for j in range(mat.shape[1]):
        #        if not mask.iloc[i, j]:
        #            val = mat.iloc[i, j]
        #            if np.abs(val) >= fc_star_thresh:
        #                ax.text(j + 0.5, i + 0.5, "*",
        #                        ha="center", va="center",
        #                        color="black", fontsize=14, fontweight="bold")

        plt.title(f"{study_col} enrichment across {cluster_label}\n"
                  f"(q < {qval_cutoff}; * = |log2FC| ≥ {fc_star_thresh})")
        plt.xlabel(study_col)
        plt.ylabel(cluster_label)
        plt.tight_layout()
        print(plt)

    return test_df, enrich_df, plt
