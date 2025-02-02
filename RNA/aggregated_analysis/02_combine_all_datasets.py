# Combine all of the fetal, postnatal non-diseased, and postnatal diseased datasets

import scanpy as sc
import numpy as np
import anndata as ad
import time

start_time = time.time()

def confirm_raw_counts(adata):
    '''
    Check if every value in adata.X.sum(axis=1), the sum of counts across all genes per cell, is an integer.
    Returns True if all sums are integers, False otherwise.
    '''
    # sum across rows (axis=1), which correspond to the genes
    sums = adata.X.sum(axis=1)

    # check if all sums are integers by comparing them to their rounded versions
    return np.all(np.equal(sums, np.round(sums)))

def convert_indices_to_string(adata):
    adata.obs.index = adata.obs.index.astype(str)
    adata.var.index = adata.var.index.astype(str)
    return adata

def process_adata(adata, metadata_to_keep):

    '''Set the adata.X as the raw counts, check that they are the raw counts, and filter adata.obs to the metadata'''

    #adata.X = adata.layers['counts']
    raw_counts_or_not = confirm_raw_counts(adata)
    print(raw_counts_or_not, flush=True)
    adata.obs = adata.obs[metadata_to_keep]

    if raw_counts_or_not:
        return(adata)
    else:
        print("Counts are not the raw data!")
        return(1)

metadata_to_keep = ["age",
                    "donor_id",
                    "sex",
                    "region",
                    "cell_type",
                    "disease",
                    "consistent_cell_type",
                    "study",
                    "technology",
                    "cell_or_nuclei", 
                    "barcode",
                    "sample_id", 
                    "age_status"]

# paths of adata files to concatenate
adata_paths = ["01_combined_LV_fetal.h5ad", "01_combined_LV_diseased.h5ad", "01_combined_LV_ND.h5ad"]

adata_list = list()

for adata_path in adata_paths:

    print(adata_path, flush=True)
    adata = sc.read_h5ad(adata_path)
    # process the adata
    adata = process_adata(adata, metadata_to_keep)
    # convert indices for all adata
    adata = convert_indices_to_string(adata)
    adata_list.append(adata)

# merge the adata together
all_adata = sc.concat(adata_list)
all_adata.layers['counts'] = all_adata.X

all_adata.write("02_combined_all_snRNA.h5ad")
print(all_adata, flush=True)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script finished in {elapsed_time}!")
