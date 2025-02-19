## Fetal reactivation spatial analysis for Kanemaru 2023

### STEP 1: Examine fetal reactivation across Visium spots

Use the script `01_analyze_spatial_fetalization_using_cell_type_proportions.ipynb` to examine the per spot, cell-type proportion weighted fetal reactivation gene signatures. 

The plots and adata files with the fetalization scores produced from this analysis are stored in `adata_with_weighted_fetalization_scores/`, created by this script. Additionally, the output `Kanemaru_Moran_I_results.csv` contains the proportion of fetalization spots and the Moran's I of fetalization.
