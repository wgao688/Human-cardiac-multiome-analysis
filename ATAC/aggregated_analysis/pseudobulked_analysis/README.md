## Perform donor pseudobulked cell-type specific analysis across 4 contrasts 

4 contrasts: sex, aging, disease, and development. Perform cell-type resolved donor-pseudobulked DAR analysis. We will only do this for the cell types that have enough (>500 nuclei) representation in the fetal, ND, and diseased groups, which include: Cardiomyocyte, Endothelial, Fibroblast, Lymphoid, Myeloid, Pericyte

### STEP 1: Perform pseudobulking
```
- 01_ATAC_peak_pseudobulker.ipynb
```

### STEP 2: Identify DARs, using a script that calls subfunctions in `run_pydeseq2_plots.py`

#### STEP 2A: Run DESeq2
```
$ nohup python3 02A_run_ATAC_pydeseq2_streamlined.py &
```

#### STEP 2B: 
Count the number of DARs per contrast 
```
- 02B_count_number_of_DARs_per_contrast.ipynb 
```

#### STEP 2C: Filter csv files to significant peaks only
As supplemental tables, filter the csv files to only the significant peaks, since otherwise this file space would be too large.

```
- 02C_filter_to_significant_DARs.ipynb
```

### STEP 3: Generate volcano plots
```
- 03_generate_volcano_plots_in_R.ipynb
```

### STEP 4: Visualize the fetalization Z-score and store the fetalization genes
```
- 04A_produce_Z_score_plot.ipynb
- 04B_examine_proportion_of_fetalization_DARs.ipynb
- 04C_plot_proportion_of_fetalization_DARs.ipynb
```

### STEP 5: Examine similarity of peaks shared across fetalization for each cell type
```
- 05_identify_fetalization_peaks.ipynb
```
