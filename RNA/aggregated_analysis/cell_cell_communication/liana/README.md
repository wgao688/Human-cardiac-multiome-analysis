## Perform liana analysis

### STEP 1: Run liana for disease vs. non-disease contrast on the scvi batch corrected adata. 

This will take about 5 hours. The input will be the final RNA adata with the scvi corrected counts. The contrast is specified by -c. This will save three outputs that will be used for the next step:

- `disease_binary_cardiac_tensor_metadata.pkl` : The metadata object 
- `disease_binary_cardiac_tensor.pkl` : The object with the actual tensor
- `01_adata_metadata.csv`: The adata.obs metadata

```
$ nohup python3 01_run_liana_for_heart_dataset_part_1.py -i ../../07_final_RNA.h5ad -c disease_binary & 
```

### STEP 2: Perform the tensor decomposition in Google Colab, so download the output files from `01_run_liana_for_heart_dataset_part_1.py` and transfer over the Google Drive. 

```
- 02_run_liana_tensor_decomposition.ipynb
```

### STEP 3: Also on Colab, perform visualization and analysis of the tensor decomposed adata.

```
- 03_analysis_post_tensor_decomposition.ipynb
```

### STEP 4: Investigate loadings for TGF-beta
```
- 04_investigation_TGF_beta_loadings.ipynb
```
