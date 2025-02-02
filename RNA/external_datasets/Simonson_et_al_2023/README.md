## Steps for preprocessing the Simonson 2023 dataset

Using single-nucleus RNA sequencing (snRNA-seq) and integrated computational analyses, Simonson et al. 2023 profiled the transcriptomes of over 99,000 human cardiac nuclei from the non-infarct region of the left ventricle of 7 ICM transplant recipients and 8 non-failing (NF) controls. 

### STEP 0: Download the dataset.

The dataset link is here: `https://singlecell.broadinstitute.org/single_cell/study/SCP1849/single-nucleus-rna-sequencing-in-ischemic-cardiomyopathy-reveals-common-transcriptional-profile-underlying-end-stage-heart-failure#/`

To downlaod the dataset with wget: 
```
wget -O ICM_scportal_05.24.2022.h5ad https://singlecell.broadinstitute.org/single_cell/data/public/SCP1849/single-nucleus-rna-sequencing-in-ischemic-cardiomyopathy-reveals-common-transcriptional-profile-underlying-end-stage-heart-failure?filename=ICM_scportal_05.24.2022.h5ad
```

### STEP 1: Preprocess the ND and D donors separately
- `01A_preprocess_Simonson_non_diseased.ipynb`
- `01B_preprocess_Simonson_diseased.ipynb`

### OUTPUTS:
- `processed_Simonson_ND.h5ad`: Processed non-diseased LV adata
- `processed_Simonson_diseased.h5ad`: Processed diseased LV adata 
