## Per-study meta-analysis for sex and aging
 
In this directory, perform analysis per study

- We will perform sex and aging analysis using pydeseq2. 
- Perform the sex analysis for datasets in which there are at least 5 female and 5 males for the given cell type.
- Perform the aging analyiss for datasets in which there are at least 5 old and 5 young for the given cell type.

### STEP 1: Run the aging analysis
```
- nohup python3 01_run_per_study_age_pydeseq2_analysis.py & 
```

this uses meta-analysis functions from `run_per_study_pydeseq2_plots.py`

### STEP 2: Run the sex analysis
```
- nohup python3 02_run_per_study_sex_pydeseq2_analysis.py &
```

### STEP 3: Examine the overlap between meta-analysis and batch-as-cov approaches

Using these scripts `03A_examine_sex_overlap.ipynb` and `03B_examine_aging_overlap.ipynb`

### STEP 4: Produce overlap plots

Using this script `04_identify_overlapping_covariate_and_meta.ipynb`
