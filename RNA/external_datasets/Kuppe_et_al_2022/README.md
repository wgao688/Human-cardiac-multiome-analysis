## Steps for preprocessing the Kuppe 2022 dataset

This dataset contains mostly myocardial infarction donors, and 4 ND donors. The healthy donors are in `00_healthy_donors.txt`. `Kuppe_supplemental_Table_1.csv` provides more information about the sampling zones for these datasets.

### STEP 0: Download dataset
Download the dataset here: https://cellxgene.cziscience.com/collections/8191c283-0816-424b-9b61-c3e1d6258a77
To get `Kuppe_et_al_snRNA.h5ad`, we can download adata version and rds with wget
```
wget -O Kuppe_et_al_snRNA.h5ad https://datasets.cellxgene.cziscience.com/c1f6034b-7973-45e1-85e7-16933d0550bc.h5ad
```

### STEP 1: Inspect the downloaded processed adata, which produces the file `01_updated_Kuppe_all.h5ad`.
- `01_inspect_all_Kuppe_adata.ipynb`

### STEP 2: Split updated adata (`01_updated_Kuppe_all.h5ad`) into ND and D, and perform preprocessing.
- `02A_preprocess_Kuppe_ND.ipynb`
- `02B_preprocess_Kuppe_D.ipynb`

### OUTPUTS:
- `02_processed_Kuppe_ND.h5ad`: Processed non-diseased LV adata
- `02_processed_Kuppe_D.h5ad`: Processed diseased LV adata
