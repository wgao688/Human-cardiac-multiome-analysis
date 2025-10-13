## Prepare inputs for scE2G

We will prepare the inputs for scE2G. It requires cell-type pseudobulked fragment files. We can either using ATAC only (~700K nuclei), 
or we can perform analysis using ~300K multiome nuclei. We will test both options. For the multiome option, we need to extract the 
RNA and ATAC for each nuclei and make them consistent. 

### STEP 1: Extract ATAC and RNA for multiome.
The interactive script `01_split_adata.ipynb` will produce paired ATAC and RNA adata objects for ENCODE and Kanemaru

### STEP 2: Confirm that barcodes are the same for each dataset, and then combine them
Perform this for Kanemaru with `02A_check_Kanemaru.ipynb`. 
Perform this for ENCODE with `02B_check_ENCODE.ipynb`.
Combine these datasets together with `02C_combine_all_ATAC_and_RNA`.

Identify cell types with at least 1 million UMIs total and 2 million fragments, based on the recommendations from scE2G

### STEP 3: Generate fragments for the different categories

- To generate for all multiome nuclei, we will perform:
```
nohup python3 03_export_fragments.py -i 03_subsampled_ATAC.h5ad -o subsampled_multiome_ATAC_fragments &
```

### STEP 4: Also export bigwig for visualization
```
nohup python3 04_export_bigwig.py -i 03_subsampled_ATAC.h5ad -o subsampled_multiome_bigwig & 
```
