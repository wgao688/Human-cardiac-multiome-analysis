# NEBULA analysis

We will test NEBULA using a fixed effect model as a strategy for identfying DEGs.

### STEP 1: Split cell types
The adata object is too large to convert as a whole to an R script, so we will run the script that sends the jobs for each cell type
```
- 01_split_adata_into_cell_types.ipynb
```

### STEP 2: Create seurat objects 

Send the `02_convert_adata_to_sce.R` script that converts adata to Seurat object for each cell type
```
nohup bash 02_send_adata_to_Seurat.sh &
```

### STEP 3: Perform NEBULA for cell-type filtered Seurat objects
```
nohup Rscript 03_iterative_NEBULA.R
```

### STEP 4: Analyze results for sex and aging
- Use fold change of 0.1 (since the fold changes are smaller for NEBULA)
Perform this interactively using `04_analyze_NEBULA_results.ipynb`
