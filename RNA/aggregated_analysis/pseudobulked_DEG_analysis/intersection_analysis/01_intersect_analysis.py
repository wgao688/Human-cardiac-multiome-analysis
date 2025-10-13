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

import intersection_scripts as funcs
from itertools import combinations

import time as time

start_time = time.time()

cell_types = ["Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial", "Epicardial",
              "Fibroblast", "Lymphoid", "LEC", "Mast", "Myeloid", "Neuronal", "Pericyte", "vSMC"]

conditions = ["aging", "disease", "development"]
contrast_pairs = list(combinations(conditions, 2))

cell_types = ["Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial", "Epicardial",
             "Fibroblast", "LEC", "Lymphoid", "Mast", "Myeloid", "Neuronal", "Pericyte", "vSMC"]

plots_dir = "intersection_plots/"
os.makedirs(plots_dir, exist_ok=True)

results_list = []  # store results here

for cell_type in cell_types:

    for contrast1, contrast2 in contrast_pairs:
        
        # load DEGs
        try: 
            contrast1_df = funcs.load_DEG_df(contrast=contrast1, cell_type=cell_type)
            contrast2_df = funcs.load_DEG_df(contrast=contrast2, cell_type=cell_type)
        
            # analyze overlaps
            up_both, down_both, up_1_down_2, down_1_up_2, consistent_ratio = funcs.analyze_gene_contrasts(
                contrast1_df=contrast1_df,
                contrast2_df=contrast2_df,
                plots_dir=plots_dir,
                cell_type=cell_type,
                contrast1=contrast1,
                contrast2=contrast2
            )

            # save overlap info
            overlap_data = {
               "cell_type": cell_type,
                "contrast1": contrast1,
                "contrast2": contrast2,
                "up_both": up_both,
                "down_both": down_both,
                "up_1_down_2": up_1_down_2,
                "down_1_up_2": down_1_up_2,
                "consistent_ratio": consistent_ratio
            }
        
            # run simulation for z-score
            z_score, p_value, z_score_plot = funcs.run_simulations(
                gene_list=list(contrast1_df.index),
                sig_df_1=contrast1_df,
                sig_df_2=contrast2_df,
                obs_ratio_of_consistent_change=consistent_ratio,
                contrast_1=contrast1,
                contrast_2=contrast2,
                num_simulations=10000,
                pseudocount=1
            )
        
            # add simulation results
            overlap_data.update({
                "z_score": z_score,
                "p_value": p_value,
                "z_score_plot": z_score_plot  # keep the plot object
            })

             # append to results list
            results_list.append(overlap_data)

        except:
             print(f"Contrast between {contrast1} and {contrast2} failed")

results_df = pd.DataFrame([{k: v for k, v in d.items() if k != "z_score_plot"} for d in results_list])
results_df.to_csv("01_intersect_analysis_results.csv")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script finished in {elapsed_time}")
