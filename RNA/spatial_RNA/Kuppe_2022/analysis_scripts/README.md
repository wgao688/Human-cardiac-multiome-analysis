## Fetal reactivation spatial analysis for Kuppe 2022

### STEP 1: Examine fetal reactivation across Visium spots

Use the script `01_analyze_spatial_fetalization_using_cell_type_proportions.ipynb` to examine the cell-type proportion weighted fetal reactivation gene signatures. 

The plots and adata files with fetalization scores are produced from this analysis are stored in `adata_with_fetalization_scores/`.

The output `Kuppe_Moran_I_results.csv` includes the proportion of fetalization spots and the Moran's I computation.

### STEP 2: Visualize cell2location cell type proportions

Use this script `02_show_cell2location_cell_type_proportions.ipynb`

### STEP 3: Model fetal reactivation score as a ordinary least squares regression against centered log transformed cell type proportions

Use this script `03_cell_type_enrichment_in_reactivation_spots.ipynb`

### STEP 4: Correlate gene expression to fetal reactivation spots

Use this script `04_gene_correlation_with_fetal_reactivation.ipynb`