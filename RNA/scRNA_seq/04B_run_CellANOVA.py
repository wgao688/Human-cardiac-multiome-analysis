import anndata as ad
import scanpy as sc
import gc
import sys
import cellanova as cnova
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sea
from collections import Counter

sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80)
pd.set_option('display.max_columns', None)
seed = 10
np.random.seed(seed)
import time

start_time = time.time()

print("Loading in the adata...", flush=True)
adata = sc.read_h5ad("03B_sc_sn_adata.h5ad")

adata_prep = cnova.model.preprocess_data(adata, integrate_key='donor_id')

# randomly sample some collection of donors across tech-plus-studyZ
adata_metadata = adata.obs
donor_level_metadata = adata_metadata[['donor_id', 'age', 'sex', 'tech_plus_study', 'cell_or_nuclei']].drop_duplicates()

sampled_donors = sampled_donors = (
    donor_level_metadata.groupby('tech_plus_study', group_keys=False)
    .apply(lambda x: x.sample(2, random_state=44))
)

print(sampled_donors.shape)

integrate_key = 'donor_id'
condition_key = 'donor_id'

control_batches = list(set(adata_prep[adata_prep.obs[condition_key].isin(sampled_donors['donor_id']),].obs[integrate_key]))

control_dict = {
    'g1': control_batches,
}
control_dict

adata_prep = cnova.model.calc_ME(adata_prep, integrate_key='donor_id')

adata_prep = cnova.model.calc_BE(adata_prep,
                                 integrate_key='donor_id',
                                 control_dict=control_dict)

adata_prep = cnova.model.calc_TE(adata_prep, integrate_key='donor_id')

adata_prep.write_h5ad('04B_post_CellANOVA.h5ad')

end_time = time.time()

elapsed_time = end_time - start_time
print(f"Elapsed time for script is {elapsed_time} s", flush=True)