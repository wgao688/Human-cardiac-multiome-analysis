import scanpy as sc 
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
import os
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from collections import Counter
from scipy.stats import rankdata

plots_dir = "plots/"
os.makedirs(plots_dir,exist_ok=True)

start_time = time.time()

# load in the adata 
adata = sc.read_h5ad("../07_final_RNA_without_scvi.h5ad")
adata

def get_filtered_adata(adata, cell_type, study=None,
                       num_cells_belonging_to_cell_type=20000, 
                       num_cells_not_belonging_to_cell_type=10000,
                      total_cells = 30000): 
    '''Get cells belonging to the input cell_type either in the consistent_cell_type or final_cell_type category. 
    To reduce compute burden later, limit cells for that cell type to a certain number using the num_cells_belonging_to_cell_type 
    parameter. For the negative control comparison, also limit it using num_cells_not_belonging_to_cell_type.'''

    # if study is specificied, then filter to just that study
    if study:
        adata = adata[adata.obs.study == study].copy()
        print(f"Filtering adata to only cells in {study}")

    cell_type_adata_in_both = adata[(adata.obs.consistent_cell_type == cell_type) & (adata.obs.final_cell_type == cell_type)]
    cell_type_adata_in_both.obs['annotation_status'] = "both"
    print(f"Cells annotated as {cell_type} in both original and revised: {cell_type_adata_in_both.shape[0]}")
    
    cell_type_only_in_original = adata[(adata.obs.consistent_cell_type == cell_type) & (adata.obs.final_cell_type != cell_type)]
    cell_type_only_in_original.obs['annotation_status'] = "original_only"
    print(f"Cells annotated as {cell_type} in only original: {cell_type_only_in_original.shape[0]}")
    
    cell_type_only_in_revised = adata[(adata.obs.consistent_cell_type != cell_type) & (adata.obs.final_cell_type == cell_type)]
    cell_type_only_in_revised.obs['annotation_status'] = "revised_only"
    print(f"Cells annotated as {cell_type} in only revised: {cell_type_only_in_revised.shape[0]}")

    # also get the cells that were annotated as neither
    cell_type_in_neither = adata[(adata.obs.consistent_cell_type != cell_type) & (adata.obs.final_cell_type != cell_type)]
    # subsample this to num_cells_not_belonging_to_cell_type cells
    cell_type_in_neither = ( cell_type_in_neither[np.random.choice(cell_type_in_neither.obs_names, 
                                                                   num_cells_not_belonging_to_cell_type, replace=False)].copy() )
    cell_type_in_neither.obs['annotation_status'] = "neither"

    # combine the adata 
    overall_adata = sc.concat([cell_type_adata_in_both, cell_type_only_in_original, cell_type_only_in_revised, cell_type_in_neither])

    overall_adata_cells = overall_adata.shape[0]

    if overall_adata_cells > total_cells:
        overall_adata = overall_adata[np.random.choice(overall_adata.obs_names, num_cells_belonging_to_cell_type, replace=False)].copy()
    
    return(overall_adata)

def calculate_pairwise_correlations(adata, count_method, correlation, annotation_key='annotation_status', hvg=False, top_genes=2000): 
    '''Calculate pairwise correlations using the log1p counts, and split the results according to annotation_key. We 
    want to get the correlations between cells belonging to a cell type. Perform this for either the scVI or raw counts.'''

    if count_method == "raw":
        # get the raw counts
        adata.X = adata.layers['counts']
        # normalize and get the log counts
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # if hvg is true, then select top_genes hvgs
        if hvg and top_genes is not None: 
            sc.pp.highly_variable_genes(adata, n_top_genes=top_genes)
            # filter to HVG
            adata = adata[:, adata.var['highly_variable']]

        # get the log normal counts
        log_norm_counts = adata.X.todense()
        
    elif count_method == "scvi":
        # get the scVI transformed counts
        log_norm_counts = adata.layers['scvi_normalized'].todense()
    else:
        print("count_method argument must be either 'raw' or 'scvi'!")
        return None

    # determine which correlation to measure
    if correlation == "pearson":
        # calculate the full correlation matrix (for Pearson)
        correlation_matrix = np.corrcoef(log_norm_counts)
    elif correlation == "spearman": 
        # rank-transform each column (Spearman correlation is based on ranks)
        ranks = np.apply_along_axis(rankdata, axis=0, arr=log_norm_counts)
        # calculate the Spearman correlation matrix
        correlation_matrix = np.corrcoef(ranks, rowvar=False)
    else: 
        print("correlation must be either 'spearman' or 'pearson'")
        return None
        
    # obtain the groups
    annotation_status = adata.obs[annotation_key].values
    groups = {
        "both": np.flatnonzero(np.array(annotation_status) == "both"),
        "revised_only": np.flatnonzero(np.array(annotation_status) == "revised_only"),
        "original_only": np.flatnonzero(np.array(annotation_status) == "original_only"),
        "neither": np.flatnonzero(np.array(annotation_status) == "neither")
    }

    # subset correlations for each group combination
    results = []
    for group1, indices1 in groups.items():
        for group2, indices2 in groups.items():
            # extract relevant correlation values
            correlations = correlation_matrix[np.ix_(indices1, indices2)].flatten()

            # remove NaNs, subsample 10,000 values if there are enough
            correlations = correlations[~np.isnan(correlations)]
            if len(correlations) > 10000:
                sampled_correlations = np.random.choice(correlations, size=10000, replace=False)
            else:
                sampled_correlations = correlations

            # store as df 
            for value in sampled_correlations:
                results.append({
                    "group1": group1,
                    "group2": group2,
                    "correlation": value
                })

    # save as df
    results_df = pd.DataFrame(results)
    
    return results_df

### make the combination order agnostic ('both:original' and 'original:both' should be the same category)
def normalize_combination(row):
    categories = sorted(row.split(':')) 
    return ':'.join(categories)  

def parse_correlation_results(results_df): 
    # get the different combinations
    results_df['combination'] = results_df['group1'] + ":" + results_df['group2']
    categories = results_df['combination'].unique()

    results_df['normalized_combination'] = results_df['combination'].apply(normalize_combination)
    results_df['normalized_combination'].unique()

    # keep only the columns of interest -- the correlations between agreed annotations, and the revised and original annotations
    # against the original annotations
    comparison_groups = ["both:both", "both:revised_only", "both:original_only", "both:neither"]
    filtered_results_df = results_df[ results_df.normalized_combination.isin(comparison_groups)]

    return(filtered_results_df)

def produce_umap_plot_original_vs_revised(filtered_adata, cell_type, save_path=None, save=False): 
    '''Produce UMAP plot of the original and revised annotations'''
    
    # filter to only those that are at least annotated as the cell type in revised or annotated
    filtered_adata =  ( filtered_adata[(filtered_adata.obs.final_cell_type == cell_type) | 
                       (filtered_adata.obs.consistent_cell_type == cell_type)] )
    sc.pl.umap(filtered_adata, color = ["final_cell_type", "consistent_cell_type"], ncols = 2, wspace=0.3, show=False)

    # make the color schemes consistent
    final_cell_types = sorted(filtered_adata.obs['final_cell_type'].unique())
    consistent_cell_types = sorted(filtered_adata.obs['consistent_cell_type'].unique())
    all_cell_types = sorted(set(final_cell_types).union(set(consistent_cell_types)))
    color_palette = sns.color_palette('tab20', len(all_cell_types))
    color_mapping = {cell_type: color_palette[i] for i, cell_type in enumerate(all_cell_types)}
    
    filtered_adata.uns['final_cell_type_colors'] = [color_mapping[cell_type] for cell_type in final_cell_types]
    filtered_adata.uns['consistent_cell_type_colors'] = [color_mapping[cell_type] for cell_type in consistent_cell_types]

    # reproduce the plt and save it 
    if save:
        with plt.rc_context():
            sc.pl.umap(filtered_adata, color=['final_cell_type', 'consistent_cell_type'], ncols = 2, wspace = 0.3, show=False)
            plt.savefig(save_path)
            plt.show()
    return plt

def produce_density_plot(filtered_results_df):
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))

    custom_colors = {
        'both:both': 'blue',
        'both:revised_only': 'black',  
        'both:original_only': 'orange',  
        'both:neither': 'gray',  
    }

    # loop through each combination and plot the KDE
    for combination, group_data in filtered_results_df.groupby('normalized_combination'):
        sns.kdeplot(
            data=group_data['correlation'], 
            fill=True,
            alpha=0.5,
            label=combination,
            color=custom_colors.get(combination, '#000000')
        )

    plt.title("Cell to cell correlations", fontsize=16)
    plt.xlabel("Pearson correlation", fontsize=14)
    plt.ylabel("density", fontsize=14)
    plt.legend(title="combination", fontsize=12, title_fontsize=14)

    return plt

def produce_violin_plot(filtered_results_df):

    sns.set(style="whitegrid")
    plt.figure(figsize=(6, 6))

    # Custom color palette
    custom_palette = {
        'both:both': 'blue',
        'both:revised_only': 'black',
        'both:original_only': 'orange',
        'both:neither': 'gray',
    }

    sns.violinplot(
        data=filtered_results_df,
        x='normalized_combination',
        y='correlation',
        palette=custom_palette,
        inner='box',  # Adds a boxplot inside the violins
        linewidth=1.2
    )
    
    plt.title("Distribution of Correlation Values by Combination", fontsize=16)
    plt.xlabel("Combination", fontsize=14)
    plt.ylabel("Pearson Correlation", fontsize=14)
    plt.xticks(rotation=45, fontsize=12)
    plt.yticks(fontsize=12)

    return plt

# run the pipeline for cell types
#cell_types = ["Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial", "Epicardial", "Fibroblast",
#             "LEC", "Lymphoid", "Mast", "Myeloid", "Neuronal", "Pericyte", "vSMC"]

rarer_cell_types = ["Adipocyte", "Endocardial", "Epicardial", "LEC", "Lymphoid", "Mast", "Neuronal", "vSMC"]
studies = ["Chaffin 2022", "Reichart 2022", "Kuppe 2022", "Kuppe 2022", "Hill 2022"]

for cell_type in rarer_cell_types:
    print(cell_type, flush=True)
    
    for study in studies: 
        print(study, flush=True)
    
        # filter full adata to that cell type and study
        filtered_adata = get_filtered_adata(adata, cell_type, study=study)
        
        # produce plot comparing the original vs. revised cell type annotations
        produce_umap_plot_original_vs_revised(filtered_adata, cell_type = cell_type, 
                                              save_path = plots_dir + cell_type + "_" + study + "_UMAP.pdf", save=True)
    
        # calculate the correlation
        results_df = calculate_pairwise_correlations(filtered_adata,
                                                 count_method = 'raw', 
                                                 correlation='pearson', 
                                                hvg=True, top_genes = 2000)
        
        # filter the results to both vs both, both vs. revised, both vs. original
        filtered_results_df = parse_correlation_results(results_df)
    
        # produce the density plot of the pairwise correlations
        density_plot = produce_density_plot(filtered_results_df)
        density_plot.savefig(plots_dir + cell_type + "_" + study + "_density_plot.pdf")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script finished in {elapsed_time} s", flush=True)
