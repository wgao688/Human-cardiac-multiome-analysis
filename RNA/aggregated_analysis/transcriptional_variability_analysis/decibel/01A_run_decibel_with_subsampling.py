import sys
import os
current_pwd = os.getcwd()
sys.path.append(current_pwd + '/module/')

import scanpy as sc
import pandas as pd
import numpy as np
import scallop as sl
import triku as tk
import matplotlib.pyplot as plt
import seaborn as sns
import math
import statsmodels.api as sm
import scanpy.external as sce
from collections import Counter
from scipy.stats import wilcoxon
from scipy.stats import ttest_ind
import decibel as dcb
import time

start_time = time.time()

# Computes metrics for all cells
# Subsamples to only n=100 cells per donor
# Performs downsampling of the counts

# load in the adata
print("Loading in the adata...", flush=True)
path_to_adata = "../../07_final_RNA_without_scvi.h5ad"
adata = sc.read_h5ad(path_to_adata)
cell_type_key = "final_cell_type"

# obtain the number of genes detected per cell and the number of transcripts detected per cell
number_of_genes_detected_per_cell = (adata.X > 0).sum(axis = 1)
number_of_transcripts_detected_per_cell = (adata.X).sum(axis = 1)
adata.obs['num_genes'] = number_of_genes_detected_per_cell
adata.obs['num_transcripts'] = number_of_transcripts_detected_per_cell

# downsample the counts
transcript_target = 500

# keep only those cells that have at least that many UMIs
adata = adata[adata.obs.num_genes > transcript_target].copy()
sc.pp.downsample_counts(adata, counts_per_cell = transcript_target)

results_dir = "decibel_results_with_subsampling/"
os.makedirs(results_dir, exist_ok = True)

# now make 'larger_cluster_cell_type' into the field called 'cell_type', for consistency with scallop
adata.obs['cell_type'] = adata.obs[cell_type_key]

# iterate through each of the cell types
cell_types = adata.obs.cell_type.unique()
min_cells_req = 100 # minimum number of cells required for analysis

for cell_type in cell_types:
    
    print(cell_type, flush=True)

    # filter the adata to that cell type
    filtered_adata = adata[adata.obs.cell_type == cell_type, :]
    
    # identify the donors and iterate
    donors = filtered_adata.obs.donor_id.unique()
    
    # keep only min_cells_req; if there are more, then subsample to min_cells_req
    indices_to_sample = []
    
    for donor in donors:
        subset_adata = filtered_adata[filtered_adata.obs.donor_id == donor, :]
        
        number_of_cells = subset_adata.obs.shape[0]
    
        if number_of_cells > min_cells_req:
            sampled_indices = np.random.choice(subset_adata.obs.index, size = min_cells_req, replace = False).tolist()
            indices_to_sample.append(sampled_indices)
        
    # flatten the list to return all of the indices to keep
    flattened_list = [value for array in indices_to_sample for value in array]

    # create a subsample adata if there are at least some cells that were sampled; otherwise skip this cell type
    if (len(flattened_list) > 0):
        sampled_adata = adata[adata.obs.index.isin(flattened_list), :]
    else:
        continue

    # perform preprocessing steps
    sc.pp.normalize_total(sampled_adata, target_sum = 1e6)
    sc.pp.log1p(sampled_adata)
    sc.pp.filter_genes(sampled_adata, min_cells = 3)

    if (sampled_adata.shape[1] <= 0):
        print("No genes left after filtering...")
        continue

    sc.pp.filter_cells(sampled_adata, min_genes = 100)
    sc.pp.highly_variable_genes(sampled_adata)
    sc.pp.pca(sampled_adata)
    sc.pp.neighbors(sampled_adata)
    tk.tl.triku(sampled_adata, use_raw=True)
    sc.pp.pca(sampled_adata)
    
    # compute the distance to centroid
    dcb.distance_to_celltype_mean(sampled_adata, batch='donor_id')
    
    # compute the distance to centroid using only cell type invariant genes
    dcb.distance_to_celltype_mean_invariant(sampled_adata, batch='donor_id')
    
    # run scallop
    scal = sl.Scallop(sampled_adata)
    sl.tl.getScore(scal, res=2.5, n_trials=30, frac_cells=0.95)
    sampled_adata.obs['freq_score'] = scal.list_bootstraps[0].freq_score
    sampled_adata.obs['scallop_noise'] = 1 - sampled_adata.obs['freq_score']
    
    # store the metadata 
    sampled_adata_metadata = sampled_adata.obs
    
    # save the file
    sampled_adata_metadata.to_csv(results_dir + cell_type + "_adata_metadata.csv", index = False)
    
    print(cell_type +  " done")

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Time elapsed for script is {elapsed_time}", flush=True)
print("Script done!")
