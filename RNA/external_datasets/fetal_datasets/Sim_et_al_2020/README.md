## Steps for preprocessing the Sim et al. 2020 datasets

This contains 9 donors, with 3 fetal, 3 young, and 3 adult donors. None of these donors were diagnosed with cardiac disease.

### STEP 0: Download data
Download the tarball from GEO with accession GSE156703. This will result in a file called `GSE156703_RAW.tar`

To unpack this tarball into individual directories for each sample, run:
```
$ nohup bash 00_process_tar.py &
```

### STEP 1: Convert to adata
We will combine the count matrices from each of the individual sample directories into a combined adata object using this interactive script.
- `01_preprocess_all_Sim.ipynb`

This will produce an intermediate adata file called `01_Sim_adata.h5ad`.

### STEP 2: Preprocess adata
We will split the adata into fetal and postnatal donors, and preprocess them separately. This is because the cutoffs for %mt, %ribo, %hb will be different for fetal and adult donors.
- `02A_preprocess_Sim_fetal.ipynb`
- `02B_preprocess_Sim_postnatal_ND.ipynb`

### OUTPUTS:
- `02_Sim_fetal.h5ad`: The adata object with fetal donors
- `02_Sim_ND.h5ad`: The adata object with non-diseased donors
