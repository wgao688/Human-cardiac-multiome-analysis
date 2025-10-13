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

def reformat_contrast_for_results_dict(contrast_name):

    '''
    Reformat a contrast name to the format for extracting it from results_dict
    Example: ("age-group", "fetal", "young") --> ['age-group_fetal_vs_young']
    '''

    factor, group1, group2 = contrast
    formatted_contrast = f"{factor}_{group1}_vs_{group2}"
    return(formatted_contrast)

cell_type = "Fibroblast"

count_matrix, metadata = funcs.load_data(cell_type = cell_type,
                                       count_matrix_dir="../raw_pseudobulked_matrices/")

contrasts = [
            ("age-group", "old", "young"),
            ("disease-binary", "Y", "N"),
            ("sex", "male", "female"),
        ]

def subsample_and_run_deseq(num_subsample, count_matrix, metadata, contrasts, funcs):
   
    # randomly sample donors
    subsampled_donors = np.random.choice(metadata.index, size=num_subsample, replace=False)

    # subset count matrix and metadata
    subsampled_metadata = metadata.loc[subsampled_donors]
    subsampled_count_matrix = count_matrix.loc[subsampled_donors]

    # run DESeq2 analysis
    results_dict, significant_genes, dds = funcs.run_deseq_analysis(
        subsampled_count_matrix,
        subsampled_metadata,
        contrasts,
        covariate_keys=["tech_plus_study", "sex", "age_group", "disease_binary"]
    )

    # create a dictionary for DEG counts
    deg_counts = {
        "num_subsample": num_subsample,
        "aging_up_DEGs": len(significant_genes['age-group_old_vs_young']['up']),
        "aging_down_DEGs": len(significant_genes['age-group_old_vs_young']['down']),
        "disease_up_DEGs": len(significant_genes['disease-binary_Y_vs_N']['up']),
        "disease_down_DEGs": len(significant_genes['disease-binary_Y_vs_N']['down']),
        "sex_up_DEGs": len(significant_genes['sex_male_vs_female']['up']),
        "sex_down_DEGs": len(significant_genes['sex_male_vs_female']['down'])
        }

    # Convert to DataFrame
    deg_df = pd.DataFrame([deg_counts])

    return deg_df

subsample_sizes = list(range(50, 290, 25)) + [290]
all_results = []

for num_subsample in subsample_sizes:
    print(num_subsample, flush=True)
    for replicate in range(5):  # 5 replicates
        print(replicate, flush=True)
        deg_results_df = subsample_and_run_deseq(num_subsample, count_matrix, metadata, contrasts, funcs)
        deg_results_df["num_subsample"] = num_subsample  # Store sample size
        deg_results_df["replicate"] = replicate + 1  # Store replicate number
        all_results.append(deg_results_df)

final_results_df = pd.concat(all_results)
final_results_df.to_csv("Fibroblast_deg_results_subsampling.csv", index=False)

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Elapsed time is {elapsed_time}")
