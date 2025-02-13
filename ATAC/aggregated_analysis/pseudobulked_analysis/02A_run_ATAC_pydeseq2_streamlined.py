import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import os
import pickle 

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data

from collections import Counter
from upsetplot import UpSet
from scipy import stats
import gseapy as gp
from gseapy import barplot, dotplot

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import run_pydeseq2_plots as funcs
import time as time

start_time = time.time()

# run for only the more common cell types that have enough cells per cell type in fetal, ND, and diseased groups
cell_types = ["Cardiomyocyte", "Endothelial", "Fibroblast", "Lymphoid", "Myeloid", "Pericyte"]

plots_dir = "DAR_contrasts_intersection_plots/"
results_dir = "pydeseq2_results/"
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

def reformat_contrast_for_results_dict(contrast_name):
    '''
    Reformat a contrast name to the format for extracting it from results_dict
    Example: ("age-group", "fetal", "young") --> ['age-group_fetal_vs_young']
    '''

    factor, group1, group2 = contrast
    formatted_contrast = f"{factor}_{group1}_vs_{group2}"
    return(formatted_contrast)

# store the Z-scores here for 3 different comparisons:
# 1. fetal & aging, 2. fetal & disease, 3. disease & aging
cell_types_fetal_disease = list()
cell_types_fetal_age = list()
cell_types_disease_age = list()

Z_scores_fetal_disease = list()
Z_scores_fetal_age = list()
Z_scores_disease_age = list()

# specify the count matrix directory
count_matrix_dir="pseudobulked_counts/"

# iterate through all of the cell types
for cell_type in cell_types:

    print(cell_type, flush=True)

    # load in the count matrix and metadata
    count_matrix, metadata = funcs.load_data(cell_type = cell_type,
                                       count_matrix_dir=count_matrix_dir)

    # identify number of fetal donors
    num_fetal_donors = metadata[metadata.age_status == "fetal"].shape[0]

    # specify the contrasts; if there are not more than 3 fetal donors, then skip that contrast
    if num_fetal_donors > 3:
        contrasts = [
            ("age-group", "fetal", "young"),
            ("age-group", "old", "young"),
            ("disease-binary", "Y", "N"),
            ("sex", "male", "female"),
        ]
    else:
         contrasts = [
            ("age-group", "old", "young"),
            ("disease-binary", "Y", "N"),
            ("sex", "male", "female"),
        ]

    # run deseq2 
    results_dict, significant_genes, dds = funcs.run_deseq_analysis(count_matrix,
                                                                    metadata, contrasts,  
                                                                    covariate_keys = ["technology", "sex", "age_group", "disease_binary"])

    # save the results_dict and dds as pickle objects
    with open(results_dir +  cell_type + "_results_dict.pkl", "wb") as f:
        pickle.dump(results_dict, f)

    with open(results_dir + cell_type + "_dds.pkl", "wb") as f:
        pickle.dump(dds, f)
    
    # produce volcano plots and save them
    for contrast in contrasts:
        formatted_contrast = reformat_contrast_for_results_dict(contrast)
    
        # produce and save volcano plot
        contrast_plot = funcs.plot_volcano(results_dict[formatted_contrast], title = formatted_contrast)
        contrast_plot.savefig(plots_dir + cell_type + "_" + formatted_contrast + "_volcano_plot.pdf")
    
        # save the results for the contrast
        contrast_results_df = results_dict[formatted_contrast]
        contrast_results_df.to_csv(results_dir + cell_type + "_" + formatted_contrast + "_results.csv")
    
    # if there are fetal samples, examine fetal-disease overlap and fetal-age overlap
    if num_fetal_donors > 3:

        # perform fetal-disease overlap
        up_both, down_both, up_1_down_2, down_1_up_2, fetal_disease_consistent_ratio = funcs.analyze_gene_contrasts(results_dict, 
                                                                                                contrast1 = 'age-group_fetal_vs_young', 
                                                                                                contrast2 = 'disease-binary_Y_vs_N', 
                                                                                                plots_dir=plots_dir, cell_type=cell_type)
        
        overlap_data = {"up_both": [up_both], "down_both": [down_both], "up_1_down_2": [up_1_down_2], "down_1_up_2": [down_1_up_2], 
                        "fetal_age_consistent_ratio": [fetal_disease_consistent_ratio]}
        overlap_df = pd.DataFrame(overlap_data)
        overlap_output_path = f"{plots_dir}/{cell_type}_fetal_disease_overlap_contrast_analysis.csv"
        overlap_df.to_csv(overlap_output_path, index=False)
    
        # obtain z-value for the significance of the overlap
        if all(var is not None for var in [up_both, down_both, up_1_down_2, down_1_up_2]) and any(var > 0 for var in [up_both, down_both, up_1_down_2, down_1_up_2]):
            z_score_fetal, p_value, z_score_plt = funcs.run_simulations(count_matrix, significance_dict = significant_genes, 
                                                 obs_ratio_of_consistent_change = fetal_disease_consistent_ratio, 
                                                 contrast_1 = 'age-group_fetal_vs_young', contrast_2 = 'disease-binary_Y_vs_N',
                                                 num_simulations = 10000, pseudocount = 1)
            z_score_plt.savefig(plots_dir + cell_type + "_disease_fetal_Z_score_simulation.pdf")
            cell_types_fetal_disease.append(cell_type)
            Z_scores_fetal_disease.append(z_score_fetal)

    
        # perform fetal-age overlap
        up_both, down_both, up_1_down_2, down_1_up_2, fetal_age_consistent_ratio = funcs.analyze_gene_contrasts(results_dict, 
                                                                                                contrast1 = 'age-group_fetal_vs_young', 
                                                                                                contrast2 = 'age-group_old_vs_young', 
                                                                                                plots_dir=plots_dir, cell_type=cell_type)

        overlap_data = {"up_both": [up_both], "down_both": [down_both], "up_1_down_2": [up_1_down_2], "down_1_up_2": [down_1_up_2],
                            "fetal_age_consistent_ratio": [fetal_age_consistent_ratio]}
        overlap_df = pd.DataFrame(overlap_data)
        overlap_output_path = f"{plots_dir}/{cell_type}_fetal_age_overlap_contrast_analysis.csv"
        overlap_df.to_csv(overlap_output_path, index=False)

        # only run if all there are DARs
        if all(var is not None for var in [up_both, down_both, up_1_down_2, down_1_up_2]) and any(var > 0 for var in [up_both, down_both, up_1_down_2, down_1_up_2]):
            # obtain z-value for the significance of the overlap
            z_score_age, p_value, z_score_plt = funcs.run_simulations(count_matrix, significance_dict = significant_genes, 
                                                 obs_ratio_of_consistent_change = fetal_age_consistent_ratio, 
                                                 contrast_1 = 'age-group_fetal_vs_young', contrast_2 = 'age-group_old_vs_young',
                                                 num_simulations = 10000, pseudocount = 1)
            z_score_plt.savefig(plots_dir + cell_type + "_fetal_age_Z_score_simulation.pdf")
            cell_types_fetal_age.append(cell_type)
            Z_scores_fetal_age.append(z_score_age)

    # also examine age + disease overlap
    up_both, down_both, up_1_down_2, down_1_up_2, age_disease_consistent_ratio = funcs.analyze_gene_contrasts(results_dict,
                                                                                                contrast1 = 'age-group_old_vs_young',
                                                                                                contrast2 = 'disease-binary_Y_vs_N',
                                                                                                plots_dir=plots_dir, 
                                                                                                cell_type=cell_type)

    overlap_data = {"up_both": [up_both], "down_both": [down_both], "up_1_down_2": [up_1_down_2], "down_1_up_2": [down_1_up_2],
                        "age_disease_consistent_ratio": [age_disease_consistent_ratio]}
    overlap_df = pd.DataFrame(overlap_data)
    overlap_output_path = f"{plots_dir}/{cell_type}_age_disease_overlap_contrast_analysis.csv"
    overlap_df.to_csv(overlap_output_path, index=False)

    if all(var is not None for var in [up_both, down_both, up_1_down_2, down_1_up_2]) and any(var > 0 for var in [up_both, down_both, up_1_down_2, down_1_up_2]):

        # obtain z-value for the significance of the overlap
        z_score_age, p_value, z_score_plt = funcs.run_simulations(count_matrix, significance_dict = significant_genes,
                                             obs_ratio_of_consistent_change = age_disease_consistent_ratio,
                                             contrast_1 = 'age-group_old_vs_young', contrast_2 = 'disease-binary_Y_vs_N',
                                             num_simulations = 10000, pseudocount = 1)

        z_score_plt.savefig(plots_dir + cell_type + "_disease_age_Z_score_simulation.pdf")
        cell_types_disease_age.append(cell_type)
        Z_scores_disease_age.append(z_score_age)

def create_z_score_plot(cell_type_list, corresponding_Z_score_list):
    '''Create z-score plot across all cell types'''

    Z_score_df = pd.DataFrame({'cell_type': cell_type_list,
                               'Z_scores': corresponding_Z_score_list})
    
    Z_score_df_sorted = Z_score_df.sort_values(by='Z_scores', ascending=False)
    
    plt.figure(figsize = (4, 6))
    sns.barplot(data = Z_score_df_sorted, x = 'Z_scores', y = 'cell_type')
    plt.axvline(1.64, color = "black", linestyle = ":")
    plt.title("Fetalization Z-scores across all cell type")
    plt.xlabel("Z-score")
    plt.ylabel("cell type")

    return(Z_score_df_sorted, plt)

print("Creating Z-score overall plots...", flush=True)
# run for fetal-age
fetal_disease_Z_df, fetal_disease_Z_fig = create_z_score_plot(cell_type_list = cell_types_fetal_disease, 
        corresponding_Z_score_list = Z_scores_fetal_disease)
fetal_disease_Z_fig.savefig(plots_dir + "Z_scores_fetal_disease.pdf", bbox_inches="tight")
fetal_disease_Z_df.to_csv(plots_dir + "Z_scores_for_fetal_disease_overlap_across_cell_types.csv")

# run for fetal-age
fetal_age_Z_df, fetal_age_Z_fig = create_z_score_plot(cell_type_list = cell_types_fetal_age, 
        corresponding_Z_score_list = Z_scores_fetal_age)
fetal_age_Z_fig.savefig(plots_dir + "Z_scores_fetal_age.pdf", bbox_inches="tight")
fetal_age_Z_df.to_csv(plots_dir + "Z_scores_for_fetal_age_overlap_across_cell_types.csv")

# run for disease-age
age_disease_Z_df, age_disease_Z_fig = create_z_score_plot(cell_type_list = cell_types_disease_age, 
        corresponding_Z_score_list = Z_scores_disease_age)
age_disease_Z_fig.savefig(plots_dir + "Z_scores_age_disease.pdf", bbox_inches="tight")
age_disease_Z_df.to_csv(plots_dir + "Z_scores_for_age_disease_overlap_across_cell_types.csv")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script finished in {elapsed_time}")
