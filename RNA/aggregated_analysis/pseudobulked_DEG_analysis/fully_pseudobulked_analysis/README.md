### Perform fully pseudobulked analysis

In this directory, perform fully pseudobulked analysis to compare to the bulk RNA-seq from GTEx results. For a fair comparison, we may want to downsample the pseudobulked counts to the same number of donors. 

### STEP 1: Run GTEx bulk RNA-seq DESeq2 for sex and aging

 - `01_run_pydeseq2_for_subsampled_GTEx.ipynb`

### STEP 2: For pseudobulked snRNA-seq, run batch as covariate approach

- `02_run_pydeseq2_batch_as_covariate.py`

### STEP 3: Por pseudobulked snRNA-seq, perform meta-analysis of studies for sex and aging
- `03_run_per_study_age_pydeseq2_analysis.py`
- `03_run_per_study_sex_pydeseq2_analysis.py`

### STEP 4: Compare results for aging and sex in terms of recall and precision.
- Use these scripts to produce and visualize results `04A_sex_compare_across_methods.ipynb` and `04B_aging_compare_across_methods.ipynb`
