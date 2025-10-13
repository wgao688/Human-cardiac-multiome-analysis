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

from upsetplot import from_contents, UpSet
import run_per_study_pydeseq2_plots as funcs
import time as time

start_time = time.time()
results_dir = "per_study_age_results_dir/"

os.makedirs(plots_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

cell_types = ["pseudobulked"]

contrasts = [
    ("age_group", "old", "young"),
]

for cell_type in cell_types:
    
    print(cell_type, flush=True)

    count_matrix, metadata = funcs.load_data(cell_type = cell_type,
                                   count_matrix_dir="../raw_pseudobulked_matrices/")

    #metadata['log_num_cells'] = np.log10(metadata['cell_count'])

    # filter to postnatal ND
    filt_metadata, filt_count_matrix = funcs.filter_to_non_diseased_postnatal(metadata = metadata, count_matrix = count_matrix)

    # identify valid studies (at least 3 in old and young)
    valid_studies = funcs.identify_valid_aging_studies(filt_metadata)

    if len(valid_studies) > 0:
        res_list = funcs.perform_deseq2_per_study(metadata = filt_metadata, count_matrix = filt_count_matrix, 
                valid_studies = valid_studies, contrasts=contrasts)

        all_age_results = funcs.collect_results(res_list, "age_group_old_vs_young")

        # run Fisher meta-analysis
        meta_results = funcs.inverse_variance_meta(df = all_age_results)

        meta_results.to_csv(results_dir + "meta_analysis_aging.csv")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script finished in {elapsed_time} s!", flush=True)
