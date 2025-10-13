## Sex and aging DEG analysis

### STEP 1: Perform DEG analysis

Perform sex and aging analysis using `01_run_pydeseq2_streamlined.py` for the "batch as covariate" approach. This calls functions in the script `run_pydeseq2_plots.py`.

### STEP 2: Extract the number of DEGs

Using this script `02_extract_sex_and_aging_DEGs.ipynb`.

### STEP 3: Produce volcano plots

Using this script `03_generate_volcano_plots_in_R.ipynb`.

### STEP 4: Compare DEG similarity across cell types.

Using these scripts `04A_extract_num_DEGs_for_covariate.ipynb` and `04B_visualize_DEG_similarity_across_cell_type.ipynb`

### STEP 5: Run gene set enrichment analysis

Using this script `05_run_GSEA_streamlined.py`. Use this script `05B_visualize_top_GSEA_R.ipynb` to visualize the results.

### Subdirectories

There are some additional subdirectories:
- `per_study_analysis`: analysis using Weighted Fisher's meta-analysis framework
- `Read_2024_analysis`: comparison to Read et al. 2024 paper
