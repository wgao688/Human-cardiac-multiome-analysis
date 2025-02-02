## Steps for preprocessing the ENCODE v4 dataset (postnatal donors)

ENCODE v4 has 54 snRNA-seq datasets which are also multiome, so they have a corresponding snATAC-seq library

### STEP 0: Download processed data 
Run this to download all of the processed adata from ENCODE server. The paths were generated using the bulk download option on the ENCODE server.

```
$ nohup xargs -L 1 curl -O -J -L < ENCODE_files.txt & 
```
This will download tar archive files. Move all of the downloaded `*.tar.gz` files after transferring into the `tarballs/` directory

### STEP 1: Extract count matrices

Each of the tar.gz files has a particular format. Enter the `tarballs/` directory and run the scripts there to process the tarballs and convert them into directories that are compatible with loading into scanpy. 

```
cd tarballs/
# extract the count matrices
$ nohup bash 01_extract_count_matrices.sh &
# identify the sample id for each of the directories by traversing the directory
$ nohup bash 02_find_ENCODE_study_number.sh & 
# prepare each of the matrices for scanpy 
$ nohup bash 03_prepare_matrix_for_scanpy.sh & 
```

Additionally, we have downloaded the metadata from ENCODE in `metadata.tsv`. Use `01_examine_and_filter_metadata.ipynb` to examine the metadata, and filter it to `filtered_metadata.tsv`.

### STEP 2: Aggregate scanpy objects together

Now that the files are in a format for scanpy, we will load them and aggregate them together into one .h5ad file using `02_create_combined_count_matrix.ipynb`. This will produce the adata file called `02_combined_ENCODE_LV_snRNA.h5ad`.

### STEP 3: Preprocess the data
We will perform standard preprocessing, as we have for the other external datasets, up to the scrublet step. This will produce the adata file called `03_prescrublet_combined_adata.h5ad`.

### STEP 4: Remove doublets with scrublet
Given the large number of nuclei, the scrublet step takes several hours. Therefore, perform it non-interactively with this python script.
```
$ nohup python3 04_run_scrublet.py
```
The resultant adata with doublets removed is `04_post_scrublet_adata.h5ad`.

### STEP 5: Examine the adata file after running scrublet interactively
- `05_examine_post_scrublet_adata.ipynb`.

### STEP 6: Annotate the cell types.
The cell types have NOT been annotated yet in the ENCODE objects. Therefore, we will perform batch integration with scVI and then annotate the cell types. This will be done in the `scVI` subdirectory: 
```
$ cd scVI
# run scVI; this will take about 2.5 hours with 99 cpus
$ nohup python3 run_scVI.py
# annotate cell types after scVI
01_visualize_and_annotate_post_scvi.ipynb
$ cd ..
```

Returning to this directory, will transfer the annotations based on the scVI batch corrected counts to `04_post_scrublet_adata.h5ad` and save it to the final adata with the cell type annotations called `post_processing_ENCODE.h5ad`. 

### OUTPUTS:
- `post_processing_ENCODE.h5ad`: The processed ENCODE LV ND adata, with cell type annotations
