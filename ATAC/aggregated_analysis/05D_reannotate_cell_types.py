import snapatac2 as snap
import scanpy as sc 
import pandas as pd
import numpy as np
import snapatac2 as snap
import matplotlib.pyplot as plt
from collections import Counter
import os
import time

start_time = time.time()

print("Loading in the adata object and metadata...", flush=True)
# load in the combined metadata
combined_metadata = pd.read_csv("05_high_quality_cell_metadata.csv", index_col = 0)
combined_metadata = combined_metadata.reset_index(drop = True)
combined_metadata.index = combined_metadata.barcode

# read in the cleaned adata
print("Reading in the original peak adata")
peak_adata = sc.read_h5ad("05_cleaned_adata.h5ad")
# make the barcode as the adata.obs index
peak_adata.obs_names = peak_adata.obs['barcode']
# ensure the order of the adata cells and the metadata are the same
peak_adata = peak_adata[combined_metadata.index , :].copy()
peak_adata.obs = peak_adata.obs.rename(columns={'barcode': 'ATAC_barcode'})

print(peak_adata)

# add cell type info 
print("Transferring cell type information...")
peak_adata.obs['cell_type'] = combined_metadata['cell_type']

# add technology and rename study to be consistent with naming in snRNA-seq
peak_adata.obs['technology'] = peak_adata.obs['study'].map(
    {
        "Penn" : "10X_ATAC",
        "ENCODE" : "Multiome-v1",
        "Kanemaru" : "Multiome-v1",
        "Kuppe" : "10X_ATAC",
    }
)

# add technology + study as an additional metadata column
peak_adata.obs['tech_plus_study'] = peak_adata.obs['technology'].astype(str) + "_" + peak_adata.obs['study'].astype(str)

peak_adata.write("05D_filtered_adata.h5ad")

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Elapsed for script in seconds: {elapsed_time}")
print("Script done!")
