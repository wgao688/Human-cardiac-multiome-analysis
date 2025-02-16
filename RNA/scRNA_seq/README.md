### Perform analysis of scRNA-seq datasets and integration of scRNA-seq with snRNA-seq

### STEP 1: Combine all scRNA-seq datasets, which include the Teichman and Koenig studies
```
- 01_combine_all_scRNA_seq.ipynb
```

### STEP 2: Run integration with scVI. This takes about 2 hours and ran for 132 epochs.
```
- nohup 02_run_scVI.py &
```

### STEP 3: Annotate the scVI adata after integration
```
- 03A_visualize_and_annotate_post_scvi.ipynb
```

### Combine the scRNA-seq with snRNA-seq data, subsampled so that there are equal number of cells and nuclei 
```
- 03B_combined_sn_sc.ipynb
```

### STEP 4: Test three different methods of integration: Harmony-only, scVI, and CellANOVA

### A: scVI
 Examine how well the single cells and nuclei overlap each other when integrated together. To do this, we will use a subsampled version of the single NUCLEUS RNA-seq dataset. Perform using either the interactive notebook or non-interactive script. 
```
- nohup python3 04A_scvi_integrate.py & (non-interactive)
- 04A_TEST_scvi.ipynb (interactive)
```

### B: CellANOVA 
```
- nohup python3 04B_run_CellANOVA.py & (non-iteractive)
- 04B_TEST_CellANOVA.ipynb 
```

# Visualize the results
```
- 04B_visualize_CellANOVA_post_integration.ipynb
```

### C: Harmony alone
```
- 04C_run_harmony_integration.ipynb
```

### D: Run no integration
```
- 04D_run_without_integration.ipynb
```

### STEP 5: Examine integration metrics with LISI
```
- 05A_examine_UMAP_embeddings.ipynb
- 05B_LISI_and_runtime_computation.ipynb
```
