# export the peak adata to bigwig format using export_coverage()

import snapatac2 as snap
import scanpy as sc
import time

start_time = time.time()

print("Load in the adata...",flush=True)
adata = sc.read_h5ad("../bigwig_bedgraph/final_adata_with_fragments.h5ad")

print("Exporting coverage to bigwigs...", flush=True)

output_dir="donor_fragment_files"

snap.ex.export_fragments(adata, groupby="sample_id", out_dir=output_dir, compression='gzip', suffix='.bed.gz')

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time for fragment file export is {elapsed_time}", flush=True)
