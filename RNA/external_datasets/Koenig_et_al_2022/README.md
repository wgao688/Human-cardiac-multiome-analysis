## Steps for preprocessing the Koenig 2022 dataset

This dataset includes ND and DCM donors. Per the paper, there are single cell n=7 (2 donor, 5 DCM) datasets, and single nuclei n=38 (25 donor, 13 DCM). This is one of the few datasets that included some scRNA-seq as well as snRNA-seq. The donor metadata is here (taken from Supplemental Table): `00_donor_metadata.txt`.

### STEP 0: Download the R object with the processed data 
GEO download link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183852

### STEP 1: Conver the Seurat obj to files that can be read into scanpy
The dataset from Koenig et al. 2022 is available in GEO as an adata object. To use this with scanpy, we need to extract the relevant information and use that to create an adata object. This script is in R, opening the object and extracting the needed info into a subdirectory called `Seurat_to_h5ad/`. From there, we will use another python script to use this information to create an adata object. 
```
Rscript 01_convert_Seurat_to_adata.R
```

### STEP 2: Convert to adata format
```
python3 02_load_into_adata.py
```

### STEP 3: Perform processing interactively to get the entire adata object with donor-level metadata
- `03_reformat_all_Koenig_adata.ipynb`

### STEP 4: Perform processing separately for diseased and non-diseased nuclei; also create the object for single cell data
- `04A_preprocess_Koenig_nuclei_non_diseased.ipynb`
- `04B_preprocess_Koenig_nuclei_diseased.ipynb`
- `04C_preprocess_Koenig_cell.ipynb`

### OUTPUTS:
- `Koenig_all_sc_snRNA.h5ad`: adata with all scRNA-seq and snRNA-seq donors
- `processed_Koenig_cell.h5ad`: adata with all scRNA-seq donors
- `processed_Koenig_diseased.h5ad`: adata with diseased snRNA-seq donors
- `processed_Koenig_ND.h5ad`: adata with ND snRNA-seq donors
