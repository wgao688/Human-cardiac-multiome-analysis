# Examines each fragment file's QC metrics

import snapatac2 as snap
import scanpy as sc
import numpy as np
import tempfile
import os
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
import shutil
# Get the temporary directory path
tmp_dir = tempfile.gettempdir()

if __name__ == '__main__':

    start_time = time.time()

    working_dir = os.getcwd() + "/"

    # load in the metadata file
    fragment_files_df = pd.read_csv("00_updated_metadata.csv", index_col = 0)
    files = fragment_files_df['full_path']
    sample_names = fragment_files_df['sample_id']

    num_samples = len(sample_names)

    temp_dir = working_dir + "tmp/"
    os.makedirs(temp_dir, exist_ok = True)

    TSSE_v_nfrag_plots_dir = "tsse_v_nfrag_plots/"
    os.makedirs(TSSE_v_nfrag_plots_dir, exist_ok = True)
    adata_dir = "prefiltering_adata_dir/"
    os.makedirs(adata_dir, exist_ok=True)

    # iterate through each of the files individually to produce the TSSE vs nfrags metrics
    for i in np.arange(num_samples):
        
        file_name = files[i]
        print(file_name, flush=True)

        sample_name = sample_names[i]
        output_h5ad = adata_dir + sample_name + ".h5ad"
    
        adata = snap.pp.import_data(fragment_file = file_name, 
                             file = output_h5ad, 
                             chrom_sizes=snap.genome.hg38,
                             sorted_by_barcode=False, 
                             tempdir = temp_dir)
        snap.metrics.tsse(adata, snap.genome.hg38)

        snap.pl.tsse(adata, out_file = TSSE_v_nfrag_plots_dir + sample_name + "_pre_trimmed_tsse_v_nfrags.png")
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Time elapsed for script is {elapsed_time}", flush = True)
