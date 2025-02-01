import scanpy as sc
import pandas as pd 
import numpy as np
import time 

start_time = time.time()
adata = sc.read_h5ad("03_prescrublet_combined_adata.h5ad")

print("Running scrublet...")

donor_key = "donor_id"

# predict the doublets with scrublet
sc.pp.scrublet(adata, batch_key = donor_key)
num_doublets = adata[adata.obs.predicted_doublet == True].shape[0]

# filter out the doublets
adata = adata[adata.obs.predicted_doublet == False, :].copy()
print(f"Filtered out {num_doublets} likely doublets")

# save adata
adata.write("04_post_scrublet_adata.h5ad")

end_time = time.time()
elapsed_time = end_time - start_time

print("Elapsed time in seconds for this script is " + str(elapsed_time) )
print("Script complete!")
