import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from upsetplot import from_contents
from upsetplot import UpSet
import pandas as pd
from upsetplot import from_contents
import gseapy as gp
import seaborn as sns
import os

results_dir = "ora_results/"
os.makedirs(results_dir, exist_ok=True)

plots_dir = "upset_plots/"
os.makedirs(plots_dir, exist_ok=True)

def identify_significant_up_and_down_DEGs(DEG_df):

    up_DEGs = DEG_df.loc[(DEG_df['significant'] == True) & (DEG_df['log2FoldChange'] > 0), :]
    down_DEGs = DEG_df.loc[(DEG_df['significant'] == True) & (DEG_df['log2FoldChange'] < 0), :]

    return([up_DEGs, down_DEGs])

def compute_overlaps(gene_dict):
    keys = list(gene_dict.keys())
    
    # only in one disease
    only_sets = {f"{k}_only": list(gene_dict[k] - set().union(*(gene_dict[kk] for kk in keys if kk != k))) for k in keys}
    
    # pairwise intersections
    pairwise_sets = {}
    for i in range(len(keys)):
        for j in range(i+1, len(keys)):
            pair_name = f"{keys[i]}_and_{keys[j]}"
            pairwise_sets[pair_name] = list((gene_dict[keys[i]] & gene_dict[keys[j]]) - set().union(*(gene_dict[kk] for kk in keys if kk not in [keys[i], keys[j]])))
    
    # shared in all three
    all_three = list(set.intersection(*gene_dict.values()))
    all_sets = {**only_sets, **pairwise_sets, "all_three": all_three}

    counts_df = pd.DataFrame({
        "Set": list(all_sets.keys()),
        "Count": [len(v) for v in all_sets.values()]
    })

    return all_sets, counts_df


def run_gsea_on_dict(overlap_dict, gene_sets=['MSigDB_Hallmark_2020'], organism="Human", cutoff=0.05):

    results_dict = {}
    all_dfs = []

    for group, genes in overlap_dict.items():
        if len(genes) == 0:
            continue  # skip empty sets

        print(f"Running enrichment for {group} ({len(genes)} genes)...")

        res = gp.enrichr(
            gene_list=genes,
            gene_sets=gene_sets,
            organism=organism,
            outdir=None,  # don’t write files
            cutoff=cutoff
        )

        # keep only significant terms
        df = res.results
        df = df[df["Adjusted P-value"] < cutoff].copy()
        df["Group"] = group

        results_dict[group] = df
        all_dfs.append(df)

    # combine all results into one dataframe
    if all_dfs:
        all_results = pd.concat(all_dfs, ignore_index=True)
    else:
        all_results = pd.DataFrame()

    return results_dict, all_results


# iterate through all cell types

cell_types = ["Adipocyte", "Cardiomyocyte", "Endothelial",
        "Fibroblast", "LEC", "Lymphoid", "Myeloid", "Neuronal", "Pericyte"]

for cell_type in cell_types:

    # Load DEG results
    DEG_files = {
        "DCM": f"pydeseq2_results/{cell_type}_disease_DCM_vs_ND_results.csv",
        "HCM": f"pydeseq2_results/{cell_type}_disease_HCM_vs_ND_results.csv",
        "ICM": f"pydeseq2_results/{cell_type}_disease_ICM_vs_ND_results.csv"
    }

    DEG_dfs = {key: pd.read_csv(file, index_col=0) for key, file in DEG_files.items()}

    # identify significant up- and down-regulated DEGs
    up_genes = {}
    down_genes = {}
    for key, df in DEG_dfs.items():
        up, down = identify_significant_up_and_down_DEGs(df)
        up_genes[key] = set(up.index)
        down_genes[key] = set(down.index)

    # compute overlaps and counts
    up_overlap_dict, up_set_sizes = compute_overlaps(up_genes)
    down_overlap_dict, down_set_sizes = compute_overlaps(down_genes)

    up_set_sizes.to_csv(results_dir + cell_type + "_upset_up_counts.csv", index=False)
    down_set_sizes.to_csv(results_dir + cell_type + "_upset_down_counts.csv", index=False)

    # upset plots
    upset_data_up = from_contents(up_genes)
    plt.figure(figsize=(8,6))
    UpSet(upset_data_up, show_counts=True, sort_by='degree').plot()
    plt.title(f"{cell_type} Upregulated DEGs")
    plt.savefig(plots_dir + cell_type + "_upset_up_DEGs.pdf")

    upset_data_down = from_contents(down_genes)
    plt.figure(figsize=(8,6))
    UpSet(upset_data_down, show_counts=True, sort_by='degree').plot()
    plt.title(f"{cell_type} Downregulated DEGs")
    plt.savefig(plots_dir + cell_type + "_upset_down_DEGs.pdf")

    # run GSEA for upregulated overlaps
    results_up_dict, all_up_results = run_gsea_on_dict(up_overlap_dict)
    all_up_results.to_csv(results_dir + cell_type + "_upset_up.csv")

    # run GSEA for downregulated overlaps
    results_down_dict, all_down_results = run_gsea_on_dict(down_overlap_dict)
    all_down_results.to_csv(results_dir + cell_type + "_upset_down.csv")
