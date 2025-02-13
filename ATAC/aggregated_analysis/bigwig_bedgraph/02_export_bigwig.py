# export the peak adata to bigwig format using export_coverage()

import snapatac2 as snap
import scanpy as sc
import time

start_time = time.time()
# get the most updated blacklist version for human genome
blacklist_path="/mnt/data1/william/human_heart_project/Final_manuscript_analysis/human_genome/hg38-blacklist.v2.bed"

print("Load in the adata...",flush=True)
adata = sc.read_h5ad("final_adata_with_fragments.h5ad")

print("Exporting coverage to bigwigs...", flush=True)

output_dir="cell_type_bigwig"

snap.ex.export_coverage(adata, 
        blacklist=blacklist_path, 
        compression=None, 
        groupby="cell_type", 
        out_dir=output_dir)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time for bigwig script is {elapsed_time}", flush=True)
