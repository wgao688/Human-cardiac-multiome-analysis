## Steps for preprocessing the Chaffin 2022 dataset

### STEP 0: Download the dataset

Download the processed data from the Single Cell Portal, requires login to this link `https://singlecell.broadinstitute.org/single_cell/study/SCP1303/single-nuclei-profiling-of-human-dilated-and-hypertrophic-cardiomyopathy#study-download`

This will result in the file with the processed adata called `human_dcm_hcm_scportal_03.17.2022.h5ad`

### STEP 1: Process the non-diseased adata interactively. This will take about 45 minutes.
```
- 01_preprocess_Chaffin_non_diseased.ipynb
```

### STEP 2: Process the diseased adata interactively. This will take about 45 minutes.
```
- 02_preprocess_Chaffin_diseased.ipynb
```

### OUTPUTS:
- `processed_Chaffin_ND.h5ad`: Processed non-diseased LV adata
- `processed_Chaffin_diseased.h5ad`: Processed diseased LV adata
