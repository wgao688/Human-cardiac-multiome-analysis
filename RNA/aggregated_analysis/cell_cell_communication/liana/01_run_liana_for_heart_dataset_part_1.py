### Given the size of the adata, run this on the server, NOT colab (which will run into RAM issues)

import pandas as pd
import scanpy as sc
import plotnine as p9
import liana as li
import cell2cell as c2c
import decoupler as dc # needed for pathway enrichment
import warnings
warnings.filterwarnings('ignore')
from collections import defaultdict
import time
import argparse
import gc

start_time = time.time()

parser = argparse.ArgumentParser(description="Run analysis on adata")
parser.add_argument('-i', '--input', type=str, required=True, help="Path to input h5ad file")
parser.add_argument('-c', '--condition_key', type=str, default='disease_binary', help="Condition key to use in analysis")

args = parser.parse_args()

# Read in adata
print("Reading in the adata..", flush=True)
adata = sc.read_h5ad(args.input)

adata_metadata = adata.obs
# save the metadata to a file 
adata_metadata = adata_metadata.to_csv("01_adata_metadata.csv")

# Filter to nuclei with enough nuclei (here set to 1000)
adata = adata[adata.obs.cell_or_nuclei == "Nuclei"]
nuclei_threshold=1000
num_nuclei_per_donor = adata.obs.groupby("donor_id").count()['study']
donors_with_enough_nuclei = num_nuclei_per_donor[num_nuclei_per_donor > nuclei_threshold].index
adata = adata[adata.obs.donor_id.isin(donors_with_enough_nuclei)].copy()

# use the scvi batch normalized counts for liana
adata.X = adata.layers['scvi_normalized']

sample_key = 'donor_id'
condition_key = args.condition_key
groupby = 'final_cell_type'

# ligand-receptor inference by sample
li.mt.rank_aggregate.by_sample(
    adata,
    groupby=groupby,
    resource_name='consensus', # NOTE: uses human gene symbols!
    sample_key=sample_key, # sample key by which we which to loop
    use_raw=False,
    verbose=True, # use 'full' to show all verbose information
    n_perms=None, # exclude permutations for speed
    return_all_lrs=True, # return all LR values
    )

sorted_samples = sorted(adata.obs[sample_key].unique())

# building a tensor
tensor = li.multi.to_tensor_c2c(adata,
                                sample_key=sample_key,
                                score_key='magnitude_rank', # can be any score from liana
                                how='outer_cells' # how to join the samples
                                )

print(tensor.tensor.shape, flush=True)

# build metadata
context_dict = adata.obs[[sample_key, condition_key]].drop_duplicates()
context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
context_dict = defaultdict(lambda: 'Unknown', context_dict)

tensor_meta = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                                  metadata_dicts=[context_dict, None, None, None],
                                                  fill_with_order_elements=True
                                                  )

# save intermediate objects
print("Saving tensor and tensor metadata...", flush=True)
c2c.io.export_variable_with_pickle(tensor, condition_key + "_cardiac_tensor.pkl")
c2c.io.export_variable_with_pickle(tensor_meta, condition_key + "_cardiac_tensor_metadata.pkl")

end_time = time.time()
elapsed_time = end_time - start_time

print(f"Elapsed time for this script is {elapsed_time}", flush=True)
print("Script complete! Proceed to step 2, which should be expedited with GPUs.", flush=True)
