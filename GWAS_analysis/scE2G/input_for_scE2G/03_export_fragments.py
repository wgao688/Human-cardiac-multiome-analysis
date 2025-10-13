# export the peak adata to bigwig format using export_coverage()

import snapatac2 as snap
import scanpy as sc
import argparse
import time

start_time = time.time()

parser = argparse.ArgumentParser(description="Export peak data to BigWig format using SnapATAC2.")
parser.add_argument('-i', '--input', required=True, help='Path to the input .h5ad file')
parser.add_argument('-o', '--output_dir', required=True, help='Path to the output directory')

args = parser.parse_args()

print("Loading the adata file...", flush=True)
output_dir = args.output_dir
adata = sc.read_h5ad(args.input)

print("Exporting coverage to bigwigs...", flush=True)

snap.ex.export_fragments(adata, 
        compression='gzip', suffix='.bed.gz', 
        groupby="cell_type", 
        out_dir=output_dir)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time for fragment script is {elapsed_time}", flush=True)
