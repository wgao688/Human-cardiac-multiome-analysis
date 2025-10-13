## Developmental DEG analysis

### STEP 1: Perform DEG analysis

Using this script `01_run_pydeseq2_development.py`, we compare fetal hearts vs. non-diseased postnatal hearts. This calls functions in this script `run_pydeseq2_plots.py`.

### STEP 2: Extract the number of DEGs

Using this script `02_extract_developmental_DEGs.ipynb`.

### STEP 3: Produce volcano plots

Using this script `03_generate_volcano_plots_in_R.ipynb`.

### STEP 4: Compare DEG similarity across cell types.

Using these scripts `04A_extract_num_DEGs_for_covariate.ipynb` and `04B_visualize_DEG_similarity_across_cell_type.ipynb`

### STEP 5: Run gene set enrichment analysis

Using this script `05_run_GSEA_streamlined.py`. Use this script `05B_visualize_top_GSEA_R.ipynb` to visualize the results.
