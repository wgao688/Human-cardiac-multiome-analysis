### Run ArchR analysis

Perform ArchR analysis on the bed files which have already been filtered to the high quality cells using SnapATAC2-based analysis, in this script's parent directory.

### STEP 1: Create metadata file for ArchR. As the bed files correspond to `sample_id`, we need to create a csv file with the metadata corresponding to each of the `sample_id` values.
```
- 01_create_ArchR_metadata.ipynb
```

### STEP 2: Create ArchR processing from the fragment files. This takes 5 hours.
```
$ nohup Rscript 02_create_ArchR_project.R &
```

### STEP 3: Perform harmony integration and actual doublet removal. This takes 5 hours.
```
$ nohup Rscript 03_further_ArchR_analysis.R
```

### STEP 4: Perform interactive visualization of UMAP and produce marker gene windows. This will produce some supplemental figures of peak-gene plots.
```
- 04_interactively_view_ArchR_project.ipynb
```

### STEP 5: Call peaks, add gene expression to ArchR, and run peak-gene correlation

# call the peaks. This takes 15 minutes.
```
- nohup Rscript 05A_call_peaks.R &
```

# add gene expression, create peak matrix, run peak2gene. This takes 15 minutes.
```
- nohup Rscript 05B_run_peak_gene_correlation.R &
```

### STEP 6: View promoter regions of marker genes
```
- 06_view_marker_genes.ipynb
```
