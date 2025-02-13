# export the peak adata to bigwig format using export_coverage()
# do this for each cell type individually, and split based on age and disease status

import snapatac2 as snap
import scanpy as sc
import time

start_time = time.time()

# get the most updated blacklist version for human genome
blacklist_path="../../../human_genome/hg38-blacklist.v2.bed"

print("Load in the adata...",flush=True)
adata = sc.read_h5ad("final_adata_with_fragments.h5ad")

print("Exporting coverage to bigwigs...", flush=True)
output_dir="cell_type_fetal_disease_bigwig"

# perform this for each cell type
cell_types = adata.obs.cell_type.unique()

for cell_type in cell_types:

    print(f"Extracting fragments for {cell_type}", flush=True)

    subset_adata = adata[adata.obs.cell_type == cell_type, :].copy()
    subset_adata.obs['disease_age_status'] = subset_adata.obs['age_status'].astype(str) + "-" + subset_adata.obs['disease_binary'].astype(str)

    # export fragments according to disease and age status
    snap.ex.export_coverage(subset_adata, 
        blacklist=blacklist_path, 
        compression=None,
        prefix = cell_type,
        groupby="disease_age_status", 
        out_dir=output_dir)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time for bigwig script is {elapsed_time}", flush=True)
