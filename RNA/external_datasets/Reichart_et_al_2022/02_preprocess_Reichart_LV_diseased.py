# Run the jupyter notebook in .py script since scrublet will take a long time

import numpy as np
import pandas as pd
import scanpy as sc
import os
from collections import Counter
import re
import gc
import scanpy.external as sce
import time

start_time = time.time() 

def convert_age_string_to_decade(age_string):
    # Extract decade number using regular expression
    
    mapping = {
        'infant stage': '0-9',
        'child stage': '0-9',
        'adolescent stage': '10-19',
        'young adult stage': '20-29',
        'third decade human stage': '30-39',
        'fourth decade human stage': '40-49',
        'fifth decade human stage': '50-59',
        'sixth decade human stage': '60-69',
        'seventh decade human stage': '70-79',
        'eighth decade human stage': '80-89'   
    }

    # Get the decade range for the provided age string
    decade = mapping.get(age_string, None)  # Returns None if age_string not found in mapping
    
    return decade

def age_convertor(input_string):
    # Split the input string by '-'
    parts = input_string.split('-')

    # Convert parts to integers
    numbers = [int(part)  for part in parts]

    # Calculate the average
    average_value = sum(numbers) / len(numbers)

    return average_value

def preprocess_adata(adata, donor_key, leiden_resolution):
    '''
    Performs: 
    1. library size normalization and log scaling
    2. identification of top 2K highly variable genes, 
    3. Principal component analysis
    4. Harmony integration
    5. Neighbor neighbors computation in Harmony integration embedding
    6. Leiden clustering 
    
    Parameters:
    adata (AnnData): adata object before preprocessing
    donor_key: the column in adata.obs that corresponds to the donor information (should be 'donor_id')
    leiden_resolution: resolution for leiden clustering, higher means more clusters will be detected

    Returns:
    adata: Postprocessed adata
    '''
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=donor_key)
    sc.tl.pca(adata)
    sce.pp.harmony_integrate(adata, donor_key)
    sc.pp.neighbors(adata, use_rep = "X_pca_harmony")
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=donor_key, size=2)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution = 0.5)
    return(adata)

# load in the full adata -- takes a while since it's quite large
adata = sc.read_h5ad("Reichart_et_al_2022.h5ad")
# filter to only those that are not designated as "normal"
adata = adata[adata.obs.disease != "normal", :]
adata = adata.raw.to_adata()
adata = adata[adata.obs.Region_x == "LV", :]
# add Reichart to make the donor id names unique
adata.obs.donor_id = adata.obs.donor_id.astype('str')
adata.obs.donor_id = "Reichart" + adata.obs.donor_id

# keep the raw counts 
adata.layers['counts'] = adata.X

# convert string from 'seventh decade human stage' to 70-79
adata.obs['updated_development_stage'] = adata.obs.development_stage.apply(convert_age_string_to_decade)

# then, convert 70-79 to 74.5 (midpoint)
adata.obs['age'] = adata.obs['updated_development_stage'].apply(lambda x: age_convertor(x))

# extract the relevant metadata
metadata = adata.obs[['cell_type', 'donor_id', 'Region_x', 'sex', 'age', 'disease']]

# for later consistency, change the column names
metadata = metadata.rename(columns={'assay': 'technology', 'Region_x': 'region'})
# add additional metadata column
metadata['technology'] = '3prime-v3'
metadata['cell_or_nuclei'] = 'Nuclei'
metadata['study'] = 'Reichart 2022'

# add back metadata to adata
adata.obs = metadata

# check that the var names are gene symbols, not ENSEMBL ID
adata.var_names
# as they are not the gene symbols, update them to the gene symbols
adata.var['ensembl_id'] = adata.var_names
adata.var_names = adata.var['feature_name']
adata.var = adata.var.drop(columns = ["feature_name"])
adata.var_names

adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hb", "pct_counts_ribo"], 
             jitter=0.4,multi_panel=True,)

# filter the cells out that are above the specified thresholds for mitochondrial, ribosomal, and hemoglobin reads
mito_threshold = 5.0
ribo_threshold = 5.0
hb_threshold = 5.0 

adata_size_before = adata.shape[0]

adata = ( adata[(adata.obs.pct_counts_mt <= mito_threshold) &
                (adata.obs.pct_counts_ribo <= ribo_threshold) &
                (adata.obs.pct_counts_hb <= hb_threshold), :].copy()
        )

adata_size_after = adata.shape[0]
num_filtered = adata_size_before - adata_size_after
print(f"Filtered out {num_filtered} cells")

donor_key = "donor_id"
# predict the doublets with scrublet
sc.pp.scrublet(adata, batch_key = donor_key)
num_doublets = adata[adata.obs.predicted_doublet == True].shape[0]
# filter out the doublets
adata = adata[adata.obs.predicted_doublet == False, :].copy()
print(f"Filtered out {num_doublets} likely doublets")

LEIDEN_RES = 0.5
adata = preprocess_adata(adata = adata, 
                         donor_key = "donor_id",
                         leiden_resolution = LEIDEN_RES)

adata.obs["consistent_cell_type"] = adata.obs["cell_type"].map(
    {
        "adipocyte": "Adipocyte",
        "cardiac muscle cell": "Cardiomyocyte",
        "cardiac neuron": "Neuronal",
        "endothelial cell": "Endothelial",
        "fibroblast of cardiac tissue": "Fibroblast",
        "lymphocyte": "Lymphoid",
        "mast cell": "Mast",
        "mural cell": "Pericyte",
        "myeloid cell": "Myeloid"
    }
)

adata = adata[adata.obs.cell_type != 'unknown', :]
adata.write("processed_Reichart_diseased.h5ad")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script finished in {elapsed_time} seconds")
print("Script complete!")
