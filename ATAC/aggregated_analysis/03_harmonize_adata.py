# Performs harmonization on the combined adata, using the sample name
# don't run the get_age_group() function since these are all fetal datasets

import snapatac2 as snap
import scanpy as sc
import numpy as np
import tempfile
import os
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time

start_time = time.time()

# Get the temporary directory path
tmp_dir = tempfile.gettempdir()

# load the combined adata file in
adata = sc.read_h5ad("02_combined_adata.h5ad")

# remove the "_QC_trimmed.h5ad_" which is still in the sample and adata.obs_names
adata.obs['sample'] = adata.obs['sample'].str.replace("_QC_trimmed.h5ad", "") 
adata.obs_names = adata.obs_names.str.replace("_QC_trimmed.h5ad", "")

# retain the original barcode information
adata.obs = adata.obs.reset_index().rename(columns = {'index': 'barcode', 'sample': 'sample_id'})

# add additional metadata
print("Adding additional metadata...", flush=True)

metadata_file = pd.read_csv("01C_QC_updated_metadata.csv", index_col = 0)

adata.obs = adata.obs.merge(metadata_file, on = "sample_id", how = "inner")
print(adata.shape)

print("Harmonizing the embedding...", flush = True)
snap.pp.harmony(adata, batch="sample_id", max_iter_harmony=20)
snap.tl.umap(adata, use_rep="X_spectral_harmony")

# repeat clustering
snap.pp.knn(adata, use_rep="X_spectral_harmony")
snap.tl.leiden(adata, resolution = 1)
snap.tl.umap(adata, use_rep="X_spectral_harmony")
snap.pl.umap(adata, color="leiden")

print("Writing the harmonized adata object to harmonized_adata.h5ad", flush=True)
adata.write("harmonized_adata.h5ad")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time for harmonization is {elapsed_time} s", flush=True)
print("Script complete!")
