# Steps for preprocessing the data from both Teichman papers (Litvinukova et al. 2020 and Kanemaru et al. 2023), which are available in Heart Atlas v2

### STEP 0: Download the datasets
The link is here: https://www.heartcellatlas.org to download the file `Global_raw_RNA_counts.h5ad`, which we will filter to LV only in the scripts below.

Altogether their data consists of 25 adult donors (14 F and 11 M), ages 20-75. See `00_patient_metadata.txt` for the relevant metadata. The exact ages are not given, so we will use the half-age value (e.g., 70-75 --> 72.5)

### STEP 1: Filter to cells or nuclei, in both cases only to the left ventricle (LV); make their cell type annotations consistent with our naming scheme
- `01A_filter_Teichman_nuclei.ipynb`
- `01B_filter_Teichman_cells.ipynb`

### STEP 2: Individually preprocess the Teichman nuclei samples that were produced using different technologies; there were three technologies: 10X v2, 10X v3, and Multiome (performed in Kanemaru 2023)
- `02A_preprocess_LV_Teichman_nuclei_GEX_v2.ipynb`
- `02B_preprocess_LV_Teichman_nuclei_GEX_v3.ipynb`
- `02C_preprocess_LV_Teichman_nuclei_Multiome.ipynb`  

### STEP 3: Recombine the nuclei adata from each technology together, post processing
- `03_combine_processed_LV_nuclei_Teichman.ipynb`

### OUTPUTS:
- `processed_10x_GEX_v2_adata.h5ad`: Processed 10X 3'v2 data
- `processed_10x_GEX_v3_adata.h5ad`: Processed 10X 3'v3 data
- `processed_Multiome_adata.h5ad`: Processed 10X multiome data
- `processed_LV_all.h5ad`: Processed all LV data
