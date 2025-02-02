# peak the ATAC adata after the peak calling has been performed
# create a file that combine the peak counts 
import scanpy as sc
import snapatac2 as snap
import time

start_time = time.time()

print("Loading in the ATAC adata after peak calling...",flush=True)
ATAC_adata = sc.read_h5ad("06_peaks.h5ad")

# merge the peaks
peaks = snap.tl.merge_peaks(ATAC_adata.uns['macs3'], snap.genome.hg38)

# create the peak matrix
peak_mat = snap.pp.make_peak_matrix(ATAC_adata, use_rep=peaks['Peaks'])

# save the peak matrix
print("Saving the peak matrix...", flush=True)
peak_mat.write("06B_peak_matrix.h5ad")

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Elapsed time for script is {elapsed_time} s",flush=True)
