import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from upsetplot import UpSet
from adjustText import adjust_text
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
from scipy import stats
from statsmodels.stats.multitest import multipletests
from scipy.stats import combine_pvalues
from upsetplot import from_contents, UpSet

def plot_volcano(results, title='Volcano Plot', fc_threshold=0.5, pval_threshold=0.05, 
                 label_threshold=0.05, max_labels=20, label_fontsize=10, eps=1e-300):
    '''
    Create a volcano plot from DESeq2 results with dynamically placed labels using adjustText.

    Inputs:
    - results: DataFrame with 'log2FoldChange', 'padj', and 'gene_id'.
    - title: Plot title.
    - fc_threshold: Fold-change threshold for significance.
    - pval_threshold: P-value threshold for significance.
    - label_threshold: P-value threshold for labeling significant genes.
    - max_labels: Maximum number of gene labels to display.
    - label_fontsize: Font size for gene labels.
    - eps: Small value to avoid log of zero.
    Returns:
    - fig: Matplotlib figure object.
    '''
    # add columns for significance and -log10(padj)
    results['-log10(padj)'] = -np.log10(results['padj'] + eps)
    results['significant'] = (abs(results['log2FoldChange']) > fc_threshold) & (results['padj'] < pval_threshold)

    # count number of significant genes
    upregulated = results[(results['log2FoldChange'] > fc_threshold) & (results['padj'] < pval_threshold)].shape[0]
    downregulated = results[(results['log2FoldChange'] < -fc_threshold) & (results['padj'] < pval_threshold)].shape[0]

    x_max = abs(results['log2FoldChange']).max()
    x_max = np.min([x_max, 10])

    # create the plot
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.scatterplot(data=results, x='log2FoldChange', y='-log10(padj)', hue='significant',
                    palette={True: 'red', False: 'black'}, alpha=0.6, ax=ax)
    ax.set_xlim(-1 * x_max, x_max)

    # add lines for the log2FC and p-value thresholds
    ax.axhline(y=-np.log10(pval_threshold), linestyle='--', color='grey')
    ax.axvline(x=fc_threshold, linestyle='--', color='grey')
    ax.axvline(x=-fc_threshold, linestyle='--', color='grey')

    # add title and labels
    ax.set_title(f"{title}\nUpregulated: {upregulated}, Downregulated: {downregulated}")
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10 Adjusted p-value')
    ax.legend(title='Significant', bbox_to_anchor=(1.05, 1), loc='upper left')

    # add dynamically placed labels for significant genes
    significant_genes = results[results['significant'] & (results['padj'] < label_threshold)]
    top_genes = significant_genes.nsmallest(max_labels, 'padj')  # Select top genes with smallest p-value

    texts = []
    for i, row in top_genes.iterrows():
        texts.append(ax.text(row['log2FoldChange'], row['-log10(padj)'], row['gene_id'], fontsize=label_fontsize))

    # use adjustText to repel labels
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

    fig.tight_layout()
    return fig

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


def filter_to_non_diseased_postnatal(metadata, count_matrix):
    '''Filter the count matrix to non-diseased postnatal donors'''

    # filter to the non-diseased postnatal
    filt_metadata = metadata.loc[(metadata.disease_binary == "N") & (metadata.age_status == "postnatal"), ]
    filt_count_matrix = count_matrix.loc[filt_metadata.index, ]

    filt_metadata.index = filt_metadata.donor_id
    filt_count_matrix.index = filt_metadata.donor_id

    return([filt_metadata, filt_count_matrix])

def identify_valid_aging_studies(metadata_df):
    '''Identify studies with at least 3 for each category specified'''

    study_age_counts = metadata_df.groupby(["study", "age_group"])["donor_id"].nunique().reset_index(name="num_donors")

    valid_studies = (study_age_counts
        .groupby("study")
        .filter(lambda df: 
                (df.loc[df["age_group"] == "old", "num_donors"].min() >= 3) and
                (df.loc[df["age_group"] == "young", "num_donors"].min() >= 3))
        ["study"]
        .unique()
    )
    return(valid_studies)


def identify_valid_sex_studies(metadata_df):
    '''Identify studies with at least 3 for each category specified'''

    study_age_counts = metadata_df.groupby(["study", "sex"])["donor_id"].nunique().reset_index(name="num_donors")

    valid_studies = (study_age_counts
        .groupby("study")
        .filter(lambda df:
                (df.loc[df["sex"] == "female", "num_donors"].min() >= 3) and
                (df.loc[df["sex"] == "male", "num_donors"].min() >= 3))
        ["study"]
        .unique()
    )
    return(valid_studies)

def perform_deseq2_per_study(metadata, count_matrix, valid_studies, contrasts):
    '''Iterates through studies to perform DESeq2 analysis'''

    res_list = list()

    for study in valid_studies:

        # filter metadata and count matrix to just that study
        study_metadata = metadata[metadata['study'] == study]
        study_matrix = count_matrix.loc[study_metadata['donor_id'], ]
        
        results_dict, significant_genes, dds = run_deseq_analysis(study_matrix,
                                                                    study_metadata,
                                                                    contrasts,
                                                                    covariate_keys = ["sex", 
                                                                                      "age_group", 
                                                                                      "log_num_cells"])
    
        res_list.append({"study": study, "results": results_dict})

    return(res_list)

def reformat_contrast_for_results_dict(contrast_name):

    '''
    Reformat a contrast name to the format for extracting it from results_dict
    Example: ("age-group", "fetal", "young") --> ['age-group_fetal_vs_young']
    '''

    factor, group1, group2 = contrast
    formatted_contrast = f"{factor}_{group1}_vs_{group2}"
    return(formatted_contrast)

def collect_results(res_list, contrast):
    '''Collect results across studies for a given contrast.'''
    dfs = []
    for entry in res_list:
        study = entry["study"]
        if contrast in entry["results"]:
            df = entry["results"][contrast].copy()
            df["study"] = study
            dfs.append(df)
    return pd.concat(dfs, axis=0)

def inverse_variance_meta(df, gene_col="gene_id"):
    results = []
    for gene, sub in df.groupby(gene_col):
        sub = sub.dropna(subset=["log2FoldChange", "lfcSE"])
        if sub.empty:
            continue

        weights = 1 / (sub["lfcSE"] ** 2)
        weighted_fc = np.sum(weights * sub["log2FoldChange"]) / np.sum(weights)

        # perform SE of combined estimate
        meta_se = np.sqrt(1 / np.sum(weights))

        # z stat and p val
        z = weighted_fc / meta_se
        p = 2 * (1 - stats.norm.cdf(abs(z)))
        results.append({
            "gene_id": gene,
            "meta_log2FC": weighted_fc,
            "meta_SE": meta_se,
            "meta_pvalue": p,
            "n_studies": len(sub)
        })

    # save as a df
    res_df = pd.DataFrame(results)

    # merge with average logFC across studies
    avg_fc = (
        df.groupby("gene_id")["log2FoldChange"]
        .mean()
        .reset_index(name="avg_log2FC")
    )

    meta_results = res_df.merge(avg_fc, on = "gene_id")

    meta_results = meta_results.sort_values(by="meta_pvalue")

    # perform BH correction
    meta_results["meta_padj"] = multipletests(meta_results["meta_pvalue"], method="fdr_bh")[1]

    return meta_results

def get_single_study_genes(all_study_df, study_name, log2FC_thresh = 0.5, padj_thresh = 0.05):
    # filter to study
    single_study_df = all_study_df[all_study_df['study'] == study_name].copy()

    up_genes = single_study_df.loc[(single_study_df['padj'] < padj_thresh) & (single_study_df['log2FoldChange'] > log2FC_thresh), 'gene_id']
    down_genes = single_study_df.loc[(single_study_df['padj'] < padj_thresh) & (single_study_df['log2FoldChange'] < -log2FC_thresh), 'gene_id']

    return(up_genes, down_genes)

def get_significant_metadata_genes(meta_df, log2FC_thresh = 0.5, padj_thresh = 0.05):
    up_genes = meta_df.loc[(meta_df['meta_padj'] < padj_thresh) & (meta_df['avg_log2FC'] > log2FC_thresh), 'gene_id']
    down_genes = meta_df.loc[(meta_df['meta_padj'] < padj_thresh) & (meta_df['avg_log2FC'] < -log2FC_thresh), 'gene_id']

    return[up_genes, down_genes]

def get_pseudobulk_sig(pseudo_df, log2FC_thresh = 0.5, padj_thresh = 0.05):

    up_genes = pseudo_df.loc[(pseudo_df['padj'] < padj_thresh) & (pseudo_df['log2FoldChange'] > log2FC_thresh), 'gene_id']
    down_genes = pseudo_df.loc[(pseudo_df['padj'] < padj_thresh) & (pseudo_df['log2FoldChange'] < -log2FC_thresh), 'gene_id']

    return(up_genes, down_genes)

def arrange_by_overlap_strength(data_df):
    data_index_df = data_df.reset_index()
    data_index_df_numeric = data_index_df.copy()

    data_index_df_numeric.loc[:, data_index_df_numeric.columns != "id"] = (
        data_index_df_numeric.drop(columns="id").astype(int)
    )

    data_index_df_numeric['col_sums'] = data_index_df_numeric.drop(columns="id").sum(axis=1)
    data_index_df_numeric.sort_values(by = "col_sums", ascending=False)

    return data_index_df_numeric
