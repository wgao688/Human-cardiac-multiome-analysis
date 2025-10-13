import snapatac2 as snap
import scanpy as sc
import os
import pandas as pd
import numpy as np
from collections import Counter
from scipy.sparse import csr_matrix
import scipy.io
import time

start_time = time.time()

# load in the count matrix
print("Reading in the SnapATAC2 adata...", flush=True)

# this peak matrix is the meta-cell aggregated peak matrix
peak_mat = sc.read_h5ad("01_subsampled_ATAC.h5ad")

peak_mat.obs.to_csv("cell_metadata.csv")
peak_mat.var.to_csv("peak_names.csv")

# transpose the dimensions for SCENIC+ compatibility
peak_mat_transpose = peak_mat.X.transpose()
print(peak_mat_transpose.shape)
peak_mat_transpose = csr_matrix(peak_mat_transpose)

# write the sparse matrix
print("Writing the sparse matrix...", flush=True)
scipy.io.mmwrite('peak_counts.mtx', peak_mat_transpose)

# also extract consensus peaks
print("Extracting consensus peaks in bed format...", flush=True)
consensus_peaks = peak_mat.var_names
consensus_peaks_df = pd.DataFrame(consensus_peaks, columns=['full_coordinates'])
# Split 'full_coordinates' into 'chr', 'start', 'end'
consensus_peaks_df[['chr', 'start_end']] = consensus_peaks_df['full_coordinates'].str.split(':', expand=True)
consensus_peaks_df[['start', 'end']] = consensus_peaks_df['start_end'].str.split('-', expand=True)
consensus_peaks_df[['start', 'end']] = consensus_peaks_df[['start', 'end']].astype(int)
consensus_peaks_df = consensus_peaks_df.drop(columns=['start_end'])
df_bed = consensus_peaks_df[['chr', 'start', 'end']].reset_index(drop = True)
df_bed.to_csv('consensus_peaks.bed', sep='\t', header=False, index=False)

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Script finished in {elapsed_time} s!")
