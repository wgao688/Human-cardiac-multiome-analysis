# Rerun UMAP and clustering using the peak matrix, rather than the tile matrix

import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import numpy as np
import snapatac2 as snap
import time

start_time = time.time()

peak_adata = sc.read_h5ad("06B_peak_matrix.h5ad")

def run_preprocessing(ATAC_adata):
    '''Perform the preprocessing steps for an ATAC adata'''

    snap.pp.select_features(ATAC_adata, n_features=50000)
    snap.tl.spectral(ATAC_adata)
    snap.pp.harmony(ATAC_adata, batch="sample_id", max_iter_harmony=20)
    snap.tl.umap(ATAC_adata, use_rep="X_spectral_harmony")
    snap.pl.umap(ATAC_adata, color="sample_id", interactive=False)

    return(ATAC_adata)

processed_peak_adata = run_preprocessing(peak_adata)

end_time = time.time()
elapsed_time = end_time - start_time

processed_peak_adata.write("06C_peak_matrix_with_new_UMAP.h5ad")

print(f"Script finished in {elapsed_time} s", flush=True)
