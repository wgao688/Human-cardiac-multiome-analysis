# calls peaks for SnapATAC2 non-interactively (as it takes a while for large data)

import snapatac2 as snap
import scanpy as sc
import numpy as np
import tempfile
import os
import pandas as pd
import argparse
import time

# function to call the peaks
def call_peaks(input_file, metadata_group, output_file_adata, temp_dir, blacklist_bed_path):

    # read the adata in
    print(f"Reading in the adata {input_file} now")
    adata = sc.read_h5ad(input_file)

    # call peaks
    print("Calling peaks with MACS3 now...")
    snap.tl.macs3(adata, groupby=metadata_group, tempdir=temp_dir, max_frag_size=200, blacklist=blacklist_bed_path, n_jobs=20)

    # merge the peaks
    #merged_peaks = snap.tl.merge_peaks(adata.uns['macs3'], chrom_sizes=snap.genome.hg38)
    
    # save the files
    print("Saving the adata to " + output_file_adata)
    adata.write(output_file_adata)
    #output_file_peaks = output_file_adata + "_merged.h5ad"
    #merged_peaks.write(output_file_peaks)

# required to support parallelism
if __name__ == '__main__':

    
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Process ATAC-seq data and call peaks.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Path to the input .h5ad file.")
    parser.add_argument('-g', '--group', type=str, required=True, help="Metadata column to group by (e.g., leiden).")
    parser.add_argument('-b', '--batch', type=str, help="Batch key (e.g. donor_id)")
    parser.add_argument('-o', '--output_file_adata', required=True, type=str, help="Path to save the output .h5ad file with peaks.")
    parser.add_argument('--temp_dir', type=str, default="/mnt/data1/william/tmp", help="Temporary directory for peak calling.")
    parser.add_argument('-r', '--blacklist', type=str, required=False, help="Path to the blacklist file")
    args = parser.parse_args()

    # run the function to call peaks
    call_peaks(args.input, args.group, args.output_file_adata, args.temp_dir, args.blacklist)

    end_time = time.time()

    elapsed_time = end_time - start_time
    print("Elapsed time in seconds is " + str(elapsed_time))
    print("Script complete!")
