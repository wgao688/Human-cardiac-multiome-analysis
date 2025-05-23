## Directory for the aggregated ATAC analysis

### STEP 0: Create consistent metadata - use the original metadata and revise the donor ids to make them consistent with snRNA-seq nomenclature
`00_original_metadata.csv` -->  `00_updated_metadata.csv`  using the interactive script `00_update_metadata.ipynb`.

### STEP 1: 

#### 1A: Perform initial visualization of QC metrics for each fragment file and determine the cutoffs. 
This takes about 6 hours. Save these prefiltered adatas in `prefiltering_adata_dir/`
```
$ nohup python3 01A_examine_QC.py &
```
This will create the QC plots in `tsse_v_nfrag_plots/`

#### 1B: Manually identify QC cutoffs.
The manually identified QC cutoffs for each sample are saved in `01B_QC_filter_per_sample.txt`. 

#### 1C: Update the metadata with QC filtering thresholds
Run  `01C_update_metadata_with_QC.ipynb`, which saves an updated metadata in `01C_QC_updated_metadata.csv`

#### 1D: Perform the QC filtering on the adata objects. 
This will create a new adata that have been QC-filtered in `post_filtering_adata_outputs/`. This takes about 4 hours (using 99 cpus).
```
$ nohup python3 01D_filter_by_QC.py &
```

### STEP 2: Combine the QC filtered adata files, perform doublet removal.
This creates an overall combined adata file called `02_combined_adata.h5ad`. This takes about 7 hours (using 99 cpus).
```
$ nohup 02_combine_filtered_adatas.py &
```

### STEP 3: Perform batch harmonization using harmony.
This will then produce a harmonized adata called `harmonized_adata.h5ad`. This takes about 2 hours.
```
$ nohup python3 03_harmonize_adata.py &
```

### STEP 4:

#### 4A: On the harmonized adata, perform gene imputation with and without MAGIC. 

MAGIC gene imputation smooths the expression out and helps with sparsity. Use them gene matrices to help with the cluster annotation. This is very time-intensive, taking 38 hours with most of this due to the MAGIC imputation steps. This step may not actually be necessary as the gene imputation doesn't perform as well as label transfer from the RNA modality. For label transferring, the 737K 10X multiome barcodes, and the correspondence between ENCODE v4 ATAC and RNA biosamples, are here: `multiome_mapping_files/`
```
$ nohup python3 04_impute_gene_expression.py &
```

#### 4B: Annotate the leiden clusters
Then use the `04B_annotate_leiden_clusters.ipynb` to perform the annotation of the leiden clusters, using label transferring from the RNA annotations and gene imputation. This takes about 1 hour to run.

### STEP 5: Add back the cell type annotations to the harmonized adata

### 5A:
Also filtered out cells belonging to unclear "Mixed" clusters. This takes about 30 minutes to run.
```
nohup python3 05A_clean_up_peak_adata.py &
```

### 5B:
Then subset to the individual cell types and perform MAGIC imputation to help confirm the cell type annotation. This will take about 30 hours. These are saved here: `pre_indiv_cell_type_adata/`
```
$ nohup python3 05B_subset_to_individual_cell_types.py &  
```
### 5C:
Interactively perform annotation and filtering of nuclei that are unlikely to belong to that cell type. This manually annotates the leiden clusters and removes cells that are likely to not belong to the cell type.
- `05C_annotate_indiv_cell_type.ipynb`

### 5D: 
Use the filtered metadata to refine the cell type annotations. This takes about 15 minutes, and produces an adata object called `05D_filtered_adata.h5ad` which includes all of the fragment files, but for which peaks have still not been called. We will call peaks using this.
```
$ nohup python3 05D_reannotate_cell_types.py &
```

### STEP 6: Call peaks on this cleaned anndata, with a blacklist file; this can take quite a while (~15 hours)

Need to use the v2 blacklist (hg38-blacklist.v2.bed), not v1 (which we tried and did not get rid of all the problematic regions that have very high ATAC coverage)!

#### 6A:
```
$ nohup python3 06_call_peaks.py -i 05D_filtered_adata.h5ad -g cell_type -o 06_peaks.h5ad -r ../../human_genome/hg38-blacklist.v2.bed &
```

#### 6B:
This saves an ATAC adata with the called peaks in -o `"06_peaks.h5ad"`. We now need to merge the peaks and create a new peak matrix adata using `06_merge_peaks.py`. This creates an adata with much less storage space required as the `adata.obsm['fragment_paired']` are no longer kept. This takes about 15 minutes.
```
$ nohup python3 06B_merge_peaks.py &
```

#### 6C:
Then, rerun harmony integration on the peak matrix to produce a final UMAP. This will take 1 hour.
```
$ nohup python3 06C_run_peak_UMAP.py & 
```

#### 6D:
Visualize the UMAP interactively and identify marker peaks for each cell type. Save these marker peaks in BED-like format here: `marker_peaks_per_cell_type` using `06D_visualize_ATAC_adata_and_produce_plots.ipynb`.

#### 6E:
Produce plots of the differential peaks and enriched motifs for each cell type, save these as pickle objects using `06E_SnapATAC2_motifs_and_marker_peaks.ipynb`

### STEP 7: Make ATAC object metadata consistent with RNA metadata. 

#### 7A:
Use `07A_make_ATAC_adata_consistent_with_RNA.ipynb` to make sure that the adata objects have consistent naming.

#### 7B:
Then, store the peaks called in bed format using `07B_get_peak_bed.ipynb`. Use ChiPSeeker to annotate the called peaks with `07C_annotate_peaks_ChipSeeker.ipynb`.

### Downstream analysis directories: enter them to see their own README.md with more information about the analysis workflow
- `metadata_plots`: additional metadata plots for main and supplemental figures + tables
- `pseudobulked_analysis`: identify pseudobulked differential peak accessibility for each cell type and contrast (age, sex, disease, development)
- `bigwig_bedgraph`: directory to create bigwig and bedgraph files for visualization
- `Spurrell_intersection_analysis`: examine the degree of intersection between the peaks called in our snATAC-seq dataset and the peaks called in Spurrell et al. 2022 cardiac enhancer paper, which perform bulk H3K27ac ChIP-seq
- `ArchR_analysis`: convert the SnapATAC2 object to a format compatible with ArchR; use ArchR for some functions such as peak-gene linkage
