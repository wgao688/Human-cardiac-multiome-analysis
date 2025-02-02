import scanpy as sc 
import pandas as pd
import numpy as np
import snapatac2 as snap
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
import os
import time

start_time=time.time()

print("Reading in the harmonized adata...", flush=True)
harmonized_adata = sc.read_h5ad("harmonized_adata.h5ad")

print("Filter to the annotated metadata...",flush=True)
labeled_metadata = pd.read_csv("04B_label_transferred_metadata.csv")

# filter to the matching barcodes
filtered_adata = harmonized_adata[harmonized_adata.obs.barcode.isin(labeled_metadata.barcode)]
filtered_adata.obs_names = filtered_adata.obs.barcode
filtered_adata = filtered_adata[labeled_metadata.barcode, :]

labeled_metadata['RNA_label_transfer_cell_type'].shape
filtered_adata.obs = filtered_adata.obs.drop(columns = "barcode")
filtered_adata.obs = filtered_adata.obs.reset_index()
filtered_adata.obs['cell_type'] = labeled_metadata['RNA_label_transfer_cell_type'].reset_index(drop = True)
# remove those in the Mixed cell type clusters
filtered_adata = filtered_adata[~filtered_adata.obs.cell_type.isin(["Mixed"])].copy()

#sc.pl.umap(filtered_adata, color = "cell_type", legend_loc = "on data")
filtered_adata.write("05_cleaned_adata.h5ad")

end_time=time.time()
elapsed_time=end_time-start_time

print(f"Script finished in {elapsed_time}", flush=True)
