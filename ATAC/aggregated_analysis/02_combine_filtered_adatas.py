# Combine the adatas that were initially filtered based on their QC metrics (TSSE and nfrags)
# Further processes these adata individually by detecting doublets and creating a tile matrix

import snapatac2 as snap
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import os 
import glob
import time

start_time = time.time()

directory = "post_filtering_adata_outputs/"

# get the full paths of the QC trimmed adata files
h5ad_files = glob.glob(os.path.join(directory, '*_QC_trimmed.h5ad'))

print(f"Number of adata files {len(h5ad_files)}")

adata_list = list()

for file_path in h5ad_files:
    print(file_path)
    adata = sc.read_h5ad(file_path)
    adata_list.append(adata)

# get the sample names by removing "_QC_trimmed.h5ad" from the file name
sample_names = [f.replace('post_filtering_adata_outputs/', '').replace('._QC_trimmed.h5ad', '') for f in h5ad_files]

print("Adding tile matrix and selecting features...", flush=True)
snap.pp.add_tile_matrix(adata_list, bin_size=5000)
snap.pp.select_features(adata_list, n_features=50000)

print("Detecting and removing doublets...", flush=True)
snap.pp.scrublet(adata_list)
snap.pp.filter_doublets(adata_list)

# create new files for the further processing adatas 
processed_directory = "post_filtering_adata_outputs/further_processed/"
os.makedirs(processed_directory, exist_ok = True)

processed_adata_paths = []
for adata, sample_name in zip(adata_list, sample_names):
    processed_adata_path = os.path.join(processed_directory, sample_name + ".h5ad")
    adata.write(processed_adata_path)
    processed_adata_paths.append(processed_adata_path)

print("Combining the adata files to an anndataset", flush=True)
data = snap.AnnDataSet(adatas = [(name, path) for name, path in zip(sample_names, processed_adata_paths)], filename = "combined_anndata_set.h5ad")

print(f'Number of cells: {data.n_obs}')
print(f'Number of unique barcodes: {np.unique(data.obs_names).size}')

# make sure that all of the barcodes are unique by adding the sample id to the barcode names
unique_cell_ids = [sa + ':' + bc for sa, bc in zip(data.obs['sample'], data.obs_names)]
data.obs_names = unique_cell_ids
assert data.n_obs == np.unique(data.obs_names).size

print("Perform dimensionality reduction", flush = True)
# Perform feature selection
snap.pp.select_features(data, n_features=50000)

# Perform spectral dimensionality reduction using Laplacian Eigenmaps
snap.tl.spectral(data)
snap.tl.umap(data)

# identify clusters in this non-batched corrected data
print("Calculating knn and performing leiden clustering", flush=True)
snap.pp.knn(data)
snap.tl.leiden(data)

# add the fragment files to data object
data.obsm['fragment_paired'] = data.adatas.obsm['fragment_paired']

# use .to_adata() to convert to .h5ad
adata_combined = data.to_adata()

# save the adata and close the anndataset to avoid data corruption
print("Writing the combined adata object", flush=True)
adata_combined.write("02_combined_adata.h5ad")
data.close()

end_time = time.time()
elapsed_time = end_time - start_time
print("The elapsed time in seconds is " + str(elapsed_time))  
print("Script complete!", flush=True)
