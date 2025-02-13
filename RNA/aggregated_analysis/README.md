# Perform combined analysis of fetal, ND, and D snRNA-seq LV across multiple external dataset and internal datasets

### STEP 1: Individually run the interactive scripts for the fetal, ND, and D, which will save a .h5ad file for each. This loads in the adata files from each of the individual studies. Preprocessing of each of these datasets individually has been done in `../external_datasets` and `../internal_datasets/`
- `01A_combine_LV_ND_adata.ipynb`
- `01B_combine_fetal_LV.ipynb`
- `01C_combine_LV_D.ipynb` 

### STEP 2: Combine all of the datasets together, non-interactively.
```
$ nohup python3 02_combine_all_datasets.py &
```

This produces a concatenated adata with all ND, D, and fetal snRNA-seq datasets called `02_combined_all_snRNA.h5ad`

### STEP 3: Inspect the metadata interactively, add additional metadata to the adata.obs, and save again as 
```
- 03_check_metadata_consistency.ipynb
```
This produces an output adata called `03_combined_all_snRNA.h5ad`

### STEP 4: Run integration with Harmony on `03_combined_all_snRNA.h5ad`; also store the result without Harmony integration. This takes about 12 hours.
```
$ nohup python3 04A_run_harmony_integration.py &
```

# visualize the clusters interactively
```
- 04B_visualize_and_annotate_post_harmony.ipynb
```

### STEP 5: Also run scVI integration and further annotation refinement.

A: First, perform scVI integration on `03_combined_all_snRNA.h5ad`. This takes about 4 hours.
```
$ nohup python3 05A_run_scVI.py & 
```

Creates a subdirectory called `scVI_model/` with the trained model

B: Manually visualize the scVI clusters, removing any leiden clusters that are primarily from just one study. This removed 3 ENCODE clusters. This takes about 1 hour.
```
This is performed in 05B_visualize_and_annotate_post_scvi.ipynb
```

C: Perform subclustering at the cell type level for the scVI adata. This takes about 3 hours
```
$ nohup python3 05C_run_subclustering.py &
```
Perform the subclustering analysis, creating a subdirectory: `scvi_pre_cleanup_subset_cell_type_adata/`

D: Manually visualize the individual cell type clusters, check marker gene expression and refine annotation
```
- 05D_filter_RNA_subclusters.ipynb
```
Filter the cells and save in `scvi_post_subclustering/`

E: Combine the refined and filtered individual cell type adata objects interactively; this will save the final RNA adata object called `05E_all_snRNA_adata.h5ad`, which includes the `scvi_normalized` layer, which stores the batch corrected counts, and a much smaller adata object called `05E_all_snRNA_adata_without_scvi.h5ad` without this layer. 

We will use this adata for integration with snATAC-seq data, and we will also use it for all downstream RNA analysis, such as differential expression analysis, cell- type proportion analysis, and cell-cell communication.
`- 05E_combine_subclustered_adata.ipynb`

### STEP 6: Compare integration between scVI and harmony. We will first visualize the UMAP embeddings. Then we will compute the LISI integation metric from Korsunsky et al. 2019 to quantify the degree of integration in terms of study (higher is better) and cell type (lower is better), as well as runtime. 
```
- 06A_visualize_UMAP_embeddings.ipynb
- 06B_LISI_and_runtime_computation.ipynb
```

### STEP 7: Finalize and subsample the adata

# A: make the metadata between the final RNA object consistent with the ATAC object. The final adata object is `07_final_RNA.h5ad`. The smaller object without scvi layer is `07_final_RNA_without_scvi.h5ad`
```
- 07A_make_RNA_metadata_consistent.ipynb
```

# B: Interactively subsample the entire adata proportionally according to cell type and donor id to about 100K cells. This subsampled adata will be useful for downstream analyses that do not require the entire adata.
```
- 07B_subsample_adata.ipynb
```

### Additional subdirectories. Please review the README.md within each of these subdirectories for more information about the analysis workflow.
- `original_annotation_vs_revised`: Analysis of the original vs. revised annotation comparison
- `subclustering_analysis`: reperform leiden clustering per cell type to obtain subclusters
- `metadata_plots`: additional metadata plots for main and supplemental figures + tables
- `pseudobulked_analysis`: identification of differentially expressed genes (DEGs) using DESeq2 and gene set enrichment analysis (GSEA)
- `cell_cell_communication`:  For cell-cell communication using the recently developed `liana tensorcell2cell`, using the scvi batch corrected counts (`07_final_RNA.h5ad`)
- `cell_type_proportion_analysis`: Cell-type proportion analysis using `scanpro / propeller` using the raw counts in `07_final_RNA_without_scvi.h5ad`.
- `senescence_analysis`: senescence analysis using multiple different gene sets, including SenMayo
- `transcriptional_variability_analysis`: transcriptional noise analysis using `scallop/decibel`
