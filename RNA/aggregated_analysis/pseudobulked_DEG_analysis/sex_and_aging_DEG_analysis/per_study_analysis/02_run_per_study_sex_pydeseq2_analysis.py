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
results_dir = "sex_results_dir/"
plots_dir = "sex_plots_dir/"

cell_types = ["Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial", "Epicardial",
              "Fibroblast", "Lymphoid", "LEC", "Mast", "Myeloid", "Neuronal", "Pericyte", "vSMC"]

# create directories for plots and results
results_dir = "sex_results_dir/"
plots_dir = "sex_plots_dir/"

os.makedirs(plots_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

contrasts = [
    ("sex", "male", "female"),
]

for cell_type in cell_types:
    
    print(cell_type, flush=True)

    count_matrix, metadata = funcs.load_data(cell_type = cell_type,
                                   count_matrix_dir="../../raw_pseudobulked_matrices/")

    metadata['log_num_cells'] = np.log10(metadata['cell_count'])

    # filter to postnatal ND
    filt_metadata, filt_count_matrix = funcs.filter_to_non_diseased_postnatal(metadata = metadata, count_matrix = count_matrix)

    # identify valid studies (at least 3 in old and young)
    valid_studies = funcs.identify_valid_sex_studies(filt_metadata)

    if len(valid_studies) > 0:
        res_list = funcs.perform_deseq2_per_study(metadata = filt_metadata, count_matrix = filt_count_matrix, 
                valid_studies = valid_studies, contrasts=contrasts)

        all_sex_results = funcs.collect_results(res_list, "sex_male_vs_female")

        # run Fisher meta-analysis
        meta_results = funcs.inverse_variance_meta(df = all_sex_results)

        pseudobulk_res = pd.read_csv("../pydeseq2_results/" + cell_type + "_sex_male_vs_female_results.csv", index_col = 0)

        # get pseudobulk significant genes
        pseudobulk_up_genes, pseudobulk_down_genes = funcs.get_pseudobulk_sig(pseudobulk_res)

        # produce upset plot
        # get significant genes from each study
        up_gene_sets = {}
        down_gene_sets = {}

        for study in valid_studies:
            up_genes, down_genes = funcs.get_single_study_genes(all_study_df=all_sex_results, study_name=study)
            up_gene_sets[study] = set(up_genes)
            down_gene_sets[study] = set(down_genes)

        # add pseudobulk
        up_gene_sets["Covariate"] = set(pseudobulk_up_genes)
        down_gene_sets["Covariate"] = set(pseudobulk_down_genes)

        # get the significant genes from meta-analysis
        meta_up_genes, meta_down_genes = funcs.get_significant_metadata_genes(meta_results)
        meta_up_genes = set(meta_up_genes)
        meta_down_genes = set(meta_down_genes)

        # prepare the dictionaries
        up_contents = {study: list(genes) for study, genes in up_gene_sets.items()}
        up_contents["Weighted Fisher"] = list(meta_up_genes)

        down_contents = {study: list(genes) for study, genes in down_gene_sets.items()}
        down_contents["Weighted Fisher"] = list(meta_down_genes)

        # convert to upset plot compatible format
        up_data = from_contents(up_contents)
        down_data = from_contents(down_contents)

        # plot up DEGs
        if not up_data.empty:
            plt.figure(figsize=(10,6))
            upset = UpSet(up_data, show_counts=True, sort_categories_by=None)
            upset.plot()
            plt.title("Shared Up DEGs")
            plt.savefig(plots_dir + cell_type + "_up_upset_plot.pdf")
            up_sig = funcs.arrange_by_overlap_strength(up_data)
            up_sig.to_csv(results_dir + cell_type + "_up_overlap_results.csv")

        # plot down DEGs
        if not down_data.empty:
            plt.figure(figsize=(10,6))
            upset = UpSet(down_data, show_counts=True, sort_categories_by=None)
            upset.plot()
            plt.title("Shared Down DEGs")
            plt.savefig(plots_dir + cell_type + "_down_upset_plot.pdf")
            down_sig = funcs.arrange_by_overlap_strength(down_data)
            down_sig.to_csv(results_dir + cell_type + "_down_overlap_results.csv")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script finished in {elapsed_time} s!", flush=True)
