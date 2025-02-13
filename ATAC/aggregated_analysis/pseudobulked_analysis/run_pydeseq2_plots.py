import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from upsetplot import UpSet
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
from scipy import stats

def plot_volcano(results, title='Volcano Plot', fc_threshold=0.5, pval_threshold=0.05, label_threshold=0.05, eps=1e-300, top_n_labels=100):
    '''
    Create a volcano plot from DESeq2 results with labels for significant genes.
    
    Inputs:
    - results: DESeq2 results as a dataframe with columns 'log2FoldChange' and 'padj'.
    - title: string for plot title
    - fc_threshold: Fold change threshold for significance.
    - pval_threshold: p-value threshold for significance.
    - label_threshold: p-value threshold for labeling significant genes.
    - eps: small value to avoid taking log of p-value of 0
    '''

    # add columns for significance and -log10(padj)
    results['-log10(padj)'] = -np.log10(results['padj'] + eps)
    results['significant'] = (abs(results['log2FoldChange']) > fc_threshold) & (results['padj'] < pval_threshold)
    
    # count number of significant genes
    upregulated = results[(results['log2FoldChange'] > fc_threshold) & (results['padj'] < pval_threshold)].shape[0]
    downregulated = results[(results['log2FoldChange'] < -fc_threshold) & (results['padj'] < pval_threshold)].shape[0]

    x_max = abs(results['log2FoldChange']).max()
    x_max = np.min([x_max, 10])
    
    # Create the plot
    plt.figure(figsize=(10, 10))
    
    sns.scatterplot(data=results, x='log2FoldChange', y='-log10(padj)', hue='significant', 
                    palette={True: 'red', False: 'black'}, alpha=0.6)
    plt.xlim(-1 * x_max, x_max)

    # add lines for the log2FC and p-value threshold
    plt.axhline(y=-np.log10(pval_threshold), linestyle='--', color='grey')
    plt.axvline(x=fc_threshold, linestyle='--', color='grey')
    plt.axvline(x=-fc_threshold, linestyle='--', color='grey')
    
    # add title
    plt.title(f"{title}\nUpregulated: {upregulated}, Downregulated: {downregulated}")
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10 Adjusted p-value')
    plt.legend(title='Significant',  bbox_to_anchor=(1.05, 1), loc='upper left')

    # identify the top N labels
    top_labels = results[results['padj'] < label_threshold].sort_values('padj').head(top_n_labels)
    
    # Add labels for significant genes
    texts = []
    for i, row in top_labels.iterrows():
        if row['significant'] and row['padj'] < label_threshold:
            texts.append(plt.text(row['log2FoldChange'], row['-log10(padj)'], row['gene_id'], fontsize=6))

    return plt

def load_data(cell_type, count_matrix_dir, TPM_threshold=1, pseudocount=1):
    '''
    Load in metadata and count matrix for cell type

    Inputs:
    cell_type: The cell type to load the data for
    count_matrix_dir: The directory for where the pseudobulked counts are located
    TPM_threshold: The threshold above which to perform DE analysis (remove the genes below this threshold)
    pseudocount: pseudocount to add, to ensure that all genes do not have 0 counts, otherwise, DESeq2 can't use geometric mean to compute library size

    Returns: 
    count_matrix: the count matrix for the specified cell type
    metadata: the corresponding metadata (donor pseudobulked)
    '''

    # load in count matrix and filter out lowly expressed genes
    count_matrix = pd.read_csv(count_matrix_dir + cell_type + "_count_matrix.csv", index_col = 0)
    
    # filter out lowly expressed genes
    counts_per_gene = count_matrix.sum(axis = 0)
    TPM_per_gene = counts_per_gene.div(counts_per_gene, axis=0) * 1e6
    genes_to_keep = count_matrix.columns[TPM_per_gene >= TPM_threshold]
    count_matrix = count_matrix[genes_to_keep]

    # add a pseudocount to avoid issue of not being able to estimate library sizes
    count_matrix = count_matrix + pseudocount

    # load in metadata and make the contrasts
    metadata = pd.read_csv(count_matrix_dir + cell_type + "_metadata.csv", index_col = 0)
    metadata = metadata.loc[count_matrix.index, :]

    return([count_matrix, metadata])

def run_deseq_analysis(count_matrix, metadata, contrasts, covariate_keys, significance_threshold=0.05, log2_fc_threshold=0.5, n_cpus=8):
    '''
    Run DESeq2 model and perform analyses across the specified contrasts.

    Inputs:
    - count_matrix: count matrix
    - metadata: metadata (used for design matrix)
    - contrasts: list of contrasts in the form [(name, level_1, level_2), ...]
    - significance_threshold: p-adjusted significance cutoff
    - covariate_keys = a list of covariates to use in the design matrix (e.g. ["sex", "disease_binary"])
    - log2_fc_threshold: log2 fold change cutoff for up/down-regulated genes
    - n_cpus: Number of CPU cores for parallel processing

    Returns:
    - results_dict: dictionary with contrast names as keys and result dataFrames as values
    - significant_genes: dictionary with contrast names as keys and up/down-regulated genes as subkeys
    '''
    
    # initialize DESeq2 model
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        counts=count_matrix,
        metadata=metadata,
        design_factors=covariate_keys,
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()

    # dictionary to store overall results and significant genes
    results_dict = {}
    significant_genes = {}
    
    # perform analyses for each contrast
    for contrast_name, level1, level2 in contrasts:
        # set up and run contrast
        stat_res = DeseqStats(dds, contrast=[contrast_name, level1, level2], inference=inference)
        stat_res.summary()
        
        # process and save results
        res_df = stat_res.results_df
        res_df['gene_id'] = res_df.index
        res_df = res_df.dropna().sort_values(by="padj")
        results_dict[f"{contrast_name}_{level1}_vs_{level2}"] = res_df
        
        # identify significant up/down-regulated genes
        up_genes = res_df[(res_df['log2FoldChange'] > log2_fc_threshold) & (res_df['padj'] < significance_threshold)].index
        down_genes = res_df[(res_df['log2FoldChange'] < -log2_fc_threshold) & (res_df['padj'] < significance_threshold)].index
        significant_genes[f"{contrast_name}_{level1}_vs_{level2}"] = {'up': up_genes, 'down': down_genes}

    return results_dict, significant_genes, dds


def run_simulations(count_matrix, significance_dict, obs_ratio_of_consistent_change, contrast_1, contrast_2,
        num_simulations = 1000, pseudocount = 1):

    '''
    Perform simulations to identify the degree of significant unidirectional overlap between up and down genes between two different 
    contrasts. This takes into account the sizes of the up and down genes per contrast, and return a p-value and Z-score via a number of simulations
    drawn from a shuffled null distribution.

    Inputs:
    - count_matrix: The count matrix, used to identify the number of genes in the "gene universe:
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

    up_in_contrast_1_size = len(significance_dict[contrast_1]['up'])
    down_in_contrast_1_size = len(significance_dict[contrast_1]['down'])

    up_in_contrast_2_size = len(significance_dict[contrast_2]['up'])
    down_in_contrast_2_size = len(significance_dict[contrast_2]['down'])

    num_genes = count_matrix.shape[1]

    # create a list to store the null distribution ratio values
    ratio_values = np.zeros(num_simulations)
    
    for i in np.arange(num_simulations):
        gene_universe = np.arange(num_genes)
        
        # simulate genes up in contrast 1 (draw from the gene universe
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

def analyze_gene_contrasts(results_dict, contrast1, contrast2, plots_dir, cell_type,
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

    # extract dataframes based on contrasts
    contrast1_df = results_dict[contrast1]
    contrast2_df = results_dict[contrast2]

    # filter genes based on fold change and adjusted p-value
    up_in_contrast1_genes = contrast1_df[(contrast1_df['log2FoldChange'] > log2fc_threshold) &
                                         (contrast1_df['padj'] < p_adj_threshold)].index

    down_in_contrast1_genes = contrast1_df[(contrast1_df['log2FoldChange'] < -log2fc_threshold) &
                                           (contrast1_df['padj'] < p_adj_threshold)].index

    up_in_contrast2_genes = contrast2_df[(contrast2_df['log2FoldChange'] > log2fc_threshold) &
                                         (contrast2_df['padj'] < p_adj_threshold)].index

    down_in_contrast2_genes = contrast2_df[(contrast2_df['log2FoldChange'] < -log2fc_threshold) &
                                           (contrast2_df['padj'] < p_adj_threshold)].index

    print(f"Up in contrast 1: {up_in_contrast1_genes.shape[0]}")
    print(f"Down in contrast 1: {down_in_contrast1_genes.shape[0]}")
    print(f"Up in contrast 2: {up_in_contrast2_genes.shape[0]}")
    print(f"Down in contrast 2: {down_in_contrast2_genes.shape[0]}")

    # create sets for upset plot
    set1 = set(up_in_contrast2_genes)
    set2 = set(down_in_contrast2_genes)
    set3 = set(up_in_contrast1_genes)
    set4 = set(down_in_contrast1_genes)
    all_genes = list(set1 | set2 | set3 | set4)
        
    # perform permutation test for intersections
    up_in_contrast1_set = set(up_in_contrast1_genes)
    down_in_contrast1_set = set(down_in_contrast1_genes)
    up_in_contrast2_set = set(up_in_contrast2_genes)
    down_in_contrast2_set = set(down_in_contrast2_genes)

    up_in_both = len(up_in_contrast2_set & up_in_contrast1_set)
    down_in_both = len(down_in_contrast2_set & down_in_contrast1_set)
    up_contrast1_down_contrast2 = len(down_in_contrast2_set & up_in_contrast1_set)
    down_contrast1_up_contrast2 = len(up_in_contrast2_set & down_in_contrast1_set)

    # perform the following only if the values are not all zero
    if all(var is not None and var > 0 for var in [up_in_both, down_in_both, up_contrast1_down_contrast2, down_contrast1_up_contrast2]):

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
    else:
        print("Number of DEGs for both pairs of contrasts were all 0, so we will not produce an upset plot")

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
