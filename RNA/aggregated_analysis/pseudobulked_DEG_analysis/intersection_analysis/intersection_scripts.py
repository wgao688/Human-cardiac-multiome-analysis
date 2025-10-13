import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import os
import pickle

from collections import Counter
from upsetplot import UpSet
from scipy import stats
import gseapy as gp
from gseapy import barplot, dotplot

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from itertools import combinations

def get_significant_genes(df, padj_threshold=0.05, lfc_threshold=0.5):
    sig = df[(df['padj'] < padj_threshold) & (df['log2FoldChange'].abs() > lfc_threshold)]
    
    up_genes = sig[sig['log2FoldChange'] > 0].index
    down_genes = sig[sig['log2FoldChange'] < 0].index
    
    return up_genes, down_genes

def load_DEG_df(contrast, cell_type):
    '''Load in the DEG for the specified contrast'''

    if contrast == "disease":
        df = ( pd.read_csv("../disease_DEG_analysis/binarized/pydeseq2_results/" + cell_type + "_disease_binary_Y_vs_N_results.csv",
                           index_col = 0) )
    elif contrast == "development":
        df = ( pd.read_csv("../developmental_DEG_analysis/pydeseq2_results/" + cell_type + "_age_group_fetal_vs_young_results.csv",
                               index_col = 0) )
    elif contrast == "sex":
        df = ( pd.read_csv("../sex_and_aging_DEG_analysis/pydeseq2_results/" + cell_type + "_sex_male_vs_female_results.csv",
                               index_col = 0) )
    elif contrast == "aging":
        df = ( pd.read_csv("../sex_and_aging_DEG_analysis/pydeseq2_results/" + cell_type + "_age_group_old_vs_young_results.csv",
                               index_col = 0) )
    return df


def analyze_gene_contrasts(contrast1_df, contrast2_df, contrast1, contrast2, plots_dir, cell_type,
                           log2fc_threshold=0.5, p_adj_threshold=0.05):

    '''
    Determine the number of genes that are directionally concordant or discordant between two different contrasts

    Inputs: 
    - results_dict: The overall DE results dictionary
    - contrast_1: The first contrast to test (e.g., 'age-group_fetal_vs_young')
    - contrast_2: The second contrast to test (e.g., 'disease-binary_Y_vs_N')
    - cell_type: cell type for the results_dict (needed for saving the upset plot)
    - plots_dir: The directory to save the plots
    - log2fc_threshold: log2FC for significance
    - p_adj_threshold: p adjusted threshold

    '''
    
    # filter genes based on fold change and adjusted p-value
    up_in_contrast1_genes, down_in_contrast1_genes = get_significant_genes(contrast1_df)
    up_in_contrast2_genes, down_in_contrast2_genes = get_significant_genes(contrast2_df)

    # create sets for upset plot
    set1 = set(up_in_contrast2_genes)
    set2 = set(down_in_contrast2_genes)
    set3 = set(up_in_contrast1_genes)
    set4 = set(down_in_contrast1_genes)
    all_genes = list(set1 | set2 | set3 | set4)

    # create sets for upset plot
    set1 = set(up_in_contrast2_genes)
    set2 = set(down_in_contrast2_genes)
    set3 = set(up_in_contrast1_genes)
    set4 = set(down_in_contrast1_genes)
    all_genes = list(set1 | set2 | set3 | set4)

    # create data frame to store the sets
    data = pd.DataFrame({
        f'{contrast2}-up': [gene in set1 for gene in all_genes],
        f'{contrast2}-down': [gene in set2 for gene in all_genes],
        f'{contrast1}-up': [gene in set3 for gene in all_genes],
        f'{contrast1}-down': [gene in set4 for gene in all_genes]
    }, index=all_genes)

    # plot the UpSet data
    upset_data = data.groupby(list(data.columns)).size()
    upset = UpSet(upset_data, sort_by='degree')
    upset.plot()
    plt.savefig(f"{plots_dir}/{cell_type}_upset_plot_{contrast1}_vs_{contrast2}.pdf")
    plt.clf()

    # perform permutation test for intersections
    up_in_contrast1_set = set(up_in_contrast1_genes)
    down_in_contrast1_set = set(down_in_contrast1_genes)
    up_in_contrast2_set = set(up_in_contrast2_genes)
    down_in_contrast2_set = set(down_in_contrast2_genes)

    up_in_both = len(up_in_contrast2_set & up_in_contrast1_set)
    down_in_both = len(down_in_contrast2_set & down_in_contrast1_set)
    up_contrast1_down_contrast2 = len(down_in_contrast2_set & up_in_contrast1_set)
    down_contrast1_up_contrast2 = len(up_in_contrast2_set & down_in_contrast1_set)

    # calculate the ratio of consistent changes to total DE genes
    total_DE_genes = up_in_both + down_in_both + up_contrast1_down_contrast2 + down_contrast1_up_contrast2
    ratio_of_consistent_change = (up_in_both + down_in_both) / total_DE_genes if total_DE_genes > 0 else 0

    # print results
    print(f"Up in both {contrast1} and {contrast2}: {up_in_both}")
    print(f"Down in both {contrast1} and {contrast2}: {down_in_both}")

    up_contrast1_down_contrast2 = len(down_in_contrast2_set & up_in_contrast1_set)
    down_contrast1_up_contrast2 = len(up_in_contrast2_set & down_in_contrast1_set)

    # calculate the ratio of consistent changes to total DE genes
    total_DE_genes = up_in_both + down_in_both + up_contrast1_down_contrast2 + down_contrast1_up_contrast2
    ratio_of_consistent_change = (up_in_both + down_in_both) / total_DE_genes if total_DE_genes > 0 else 0

    # print results
    print(f"Up in both {contrast1} and {contrast2}: {up_in_both}")
    print(f"Down in both {contrast1} and {contrast2}: {down_in_both}")
    print(f"Up in {contrast1}, down in {contrast2}: {up_contrast1_down_contrast2}")
    print(f"Down in {contrast1}, up in {contrast2}: {down_contrast1_up_contrast2}")
    print("Ratio of consistent change:", ratio_of_consistent_change)

    return [up_in_both, down_in_both, up_contrast1_down_contrast2, down_contrast1_up_contrast2, ratio_of_consistent_change]

def run_simulations(gene_list, sig_df_1, sig_df_2, obs_ratio_of_consistent_change, contrast_1, contrast_2,
        num_simulations = 1000, pseudocount = 1):

    '''
    Perform simulations to identify the degree of significant unidirectional overlap between up and down genes between two different
    contrasts. This takes into account the sizes of the up and down genes per contrast, and return a p-value and Z-score via a number of simulations
    drawn from a shuffled null distribution.

    Inputs:
    - list of all genes (in the gene universe)
    - significance_dict: The dictionary with the significant results
    - ratio_of_significance_change: The actual observed degree of significant change between the two contrasts
    - contrast_1: The name of the first contrast (e.g., 'group_fetal_vs_young')
    - contrast_2: The name of the second contrast (e.g., 'disease-binary_Y_vs_N')
    - pseudocount: Value to add to the denominator of the observed to avoid dividing by 0

    Returns:
    - z_score: The z-score for the observed ratio of overlap
    - p_value: The p-values for the observed ratio of overlap
    - plt: The plot of the null distribution and the observed ratio of overlap
    '''

    sig_1_up, sig_1_down = get_significant_genes(sig_df_1)
    up_in_contrast_1_size = len(sig_1_up)
    down_in_contrast_1_size = len(sig_1_down)

    sig_2_up, sig_2_down = get_significant_genes(sig_df_2)
    up_in_contrast_2_size = len(sig_2_up)
    down_in_contrast_2_size = len(sig_2_down)

    num_genes = len(gene_list)

    # create a list to store the null distribution ratio values
    ratio_values = np.zeros(num_simulations)

    for i in np.arange(num_simulations):
        gene_universe = np.arange(num_genes)

        # simulate genes up in contrast 1 (draw from the gene universe)
        up_in_contrast1_sim = np.random.choice(gene_universe, size=up_in_contrast_1_size, replace=False)

        # simulate genes down in contrast 1: draw from remaining genes
        remaining_indices = np.setdiff1d(gene_universe, up_in_contrast1_sim)
        down_in_contrast1_sim = np.random.choice(remaining_indices, size=down_in_contrast_1_size, replace=False)

        # simulate genes up in contrast 2: draw from the whole gene universe
        up_in_contrast2_sim = np.random.choice(gene_universe, size=up_in_contrast_2_size, replace=False)

        # simulate genes down in contrast 2: draw from remaining genes
        remaining_indices = np.setdiff1d(gene_universe, down_in_contrast_2_size)
        down_in_contrast2_sim = np.random.choice(remaining_indices, size=down_in_contrast_2_size, replace=False)

        # calculate the number up in both
        n_gene_unidirectional = len(set(up_in_contrast1_sim) & set(up_in_contrast2_sim)) + len(set(down_in_contrast1_sim) & set(down_in_contrast2_sim))
        n_gene_non_directional = len(set(up_in_contrast1_sim) & set(down_in_contrast2_sim)) + len(set(down_in_contrast1_sim) & set(up_in_contrast2_sim))

        ratio_unidirectional = (n_gene_unidirectional) / (n_gene_unidirectional + n_gene_non_directional + pseudocount)
        ratio_values[i] = ratio_unidirectional

    # calculate the mean for the null distribution
    mean = np.mean(ratio_values)
    # calculate the standard deviation for the null distribution
    std_dev = np.std(ratio_values)

    # calculate the Z-score for the observed ratio of consistent change against the null distribution
    v = obs_ratio_of_consistent_change
    z_score = (v - mean) / std_dev

    # calculate the p-value for the observed value, for a two-tailed test
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

    # print the results
    print(f"Mean (μ): {mean}")
    print(f"Standard Deviation (σ): {std_dev}")
    print(f"Observed ratio of consistent change between the two contrasts (v): {v}")
    print(f"Z-score: {z_score}")
    print(f"P-value: {p_value}")

    # specify a floor on the p-value
    if p_value < 1e-16:
        p_value = 1e-16

    plt.hist(ratio_values, bins = 100)
    plt.axvline(obs_ratio_of_consistent_change,
        color='black', linestyle='--',
        label=f'Ratio of Consistent Change = {obs_ratio_of_consistent_change:.2f}')

    x_min, x_max = plt.xlim()

    # annotate the plot with the Z-score and p-value in the top-right corner
    plt.text(x_max - 0.1 * (x_max - x_min), plt.ylim()[1] * 0.95,
     f'Z-score: {z_score:.2f}\nP-value: {p_value:.2e}',
     horizontalalignment='right', verticalalignment='top',
     bbox=dict(facecolor='white', alpha=0.5))

    plt.title(f"Proportion of consistent {contrast_1} + {contrast_2} DE genes / total")
    plt.xlabel("(unidirectional DEGs) / \n (union of all DEGs)")
    plt.ylabel("number of permutations")
    plt.show()

    return([z_score, p_value, plt])
