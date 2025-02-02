# Performs QC preprocessing and initial clustering before batch correction
# Performs this using the .h5ad in memory to expedite the process

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

    filtered_adata_output_dir = "post_filtering_adata_outputs/"
    if os.path.exists(filtered_adata_output_dir):
        shutil.rmtree(filtered_adata_output_dir)
    os.makedirs(filtered_adata_output_dir)

    # load in the metadata file
    fragment_files_df = pd.read_csv("01C_QC_updated_metadata.csv", index_col = 0)

    # extract the full path
    files = fragment_files_df['full_path'].tolist()

    # extract the sample_ids
    sample_names = fragment_files_df['sample_id'].tolist()

    # use the file names to determine the output h5ad file path
    output_h5ads = [filtered_adata_output_dir + "/" + sample_name + ".h5ad" for sample_name in sample_names]
    temp_dir = working_dir + "tmp/"
    os.makedirs(temp_dir, exist_ok = True)

    adatas = snap.pp.import_data(fragment_file = files, file = output_h5ads, 
                             chrom_sizes=snap.genome.hg38,
                             sorted_by_barcode=False,
                             tempdir = temp_dir)
    snap.metrics.tsse(adatas, snap.genome.hg38)

    TSSE_v_nfrag_plots_dir = "post_filtering_tsse_v_nfrag_plots/"
    os.makedirs(TSSE_v_nfrag_plots_dir, exist_ok = True)

    num_samples = len(adatas)
    print("Perform QC filtering using predetermined filters from part 1", flush=True)

    # nfrags and tsse cutoffs
    n_frag_filters = fragment_files_df['nfrag'].tolist()
    tsse_filters = fragment_files_df['tsse'].tolist()

    adata_filtered_list = list()

    # for each of the samples, perform the QC filtering
    for i in np.arange(num_samples):
    
        sample_name = sample_names[i]
        n_frag_filter = n_frag_filters[i]
        tsse_filter = tsse_filters[i]

        print(sample_name, flush = True)

        # if there are specified cutoffs (not NA)
        print(sample_name, flush = True)
        adata_filtered = adatas[i].to_memory()
        snap.pp.filter_cells(adata_filtered, min_counts = n_frag_filter, min_tsse = tsse_filter)
        
        # display and save the plot after QC filtering
        snap.pl.tsse(adata_filtered, out_file = TSSE_v_nfrag_plots_dir + sample_names[i] + "_post_trimmed_tsse_v_nfrags.png")
        # append to the growing list of adata files
        adata_filtered_list.append(adata_filtered)
        adata_filtered.write(filtered_adata_output_dir + sample_names[i] + "_QC_trimmed.h5ad")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Time elapsed for script is " + str(elapsed_time), flush = True)
    print("script done!", flush = True)
