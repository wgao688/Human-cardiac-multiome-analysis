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

# create directories for plots and results
plots_dir = "batch_as_covariate_plots/"
results_dir = "batch_as_covariate_pydeseq2_results/"

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

cell_type = "pseudobulked"
print(cell_type, flush=True)

# load in the count matrix and metadata
count_matrix, metadata = funcs.load_data(cell_type = cell_type,
                                       count_matrix_dir="../raw_pseudobulked_matrices/")

#metadata['log_num_cells'] = np.log10(metadata['cell_count'])
    
# filter to the non-diseased postnatal
filt_metadata = metadata.loc[(metadata.disease_binary == "N") & (metadata.age_status == "postnatal"), ]
filt_count_matrix = count_matrix.loc[filt_metadata.index, ]

contrasts = [
    ("age_group", "old", "young"),
    ("sex", "male", "female"),
]

# run deseq2 
results_dict, significant_genes, dds = funcs.run_deseq_analysis(filt_count_matrix,
                                                                    filt_metadata, 
                                                                    contrasts,  
                                                                    #covariate_keys = ["tech_plus_study", "sex", "age_group", "log_num_cells"]
                                                                    covariate_keys = ["tech_plus_study", "sex", "age_group"])

# save the results_dict and dds as pickle objects
with open(results_dir +  cell_type + "_results_dict.pkl", "wb") as f:
    pickle.dump(results_dict, f)

with open(results_dir + cell_type + "_dds.pkl", "wb") as f:
    pickle.dump(dds, f)
    
# produce volcano plots and save them
for contrast in contrasts:
    print(f"Creating volcano plot for {cell_type} and contrast {contrast}", flush=True)
    formatted_contrast = reformat_contrast_for_results_dict(contrast)
    
    # produce and save volcano plot
    contrast_plot = funcs.plot_volcano(results_dict[formatted_contrast], title = formatted_contrast)
    contrast_plot.savefig(plots_dir + cell_type + "_" + formatted_contrast + "_volcano_plot.pdf")
    plt.close(contrast_plot)
    
    # save the results for the contrast
    contrast_results_df = results_dict[formatted_contrast]
    contrast_results_df.to_csv(results_dir + cell_type + "_" + formatted_contrast + "_results.csv")
    
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script finished in {elapsed_time}")
