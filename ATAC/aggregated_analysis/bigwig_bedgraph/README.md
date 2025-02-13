## Export fragments to bigwig format

This directory involves scripts to export the fragment file information to bigwig format for visualization using IGV. 

### STEP 1: Transfer the `fragments_paired` stored in the adata back to the final adata object

We need to transfer back the insertion information using this script.
`- 01_add_peaks_to_final_adata.ipynb`

### STEP 2: Export the reads in bigwig format for each cell type. This takes about 2.5 hours
```
$ nohup python3 02_export_bigwig.py &
```
### STEP 3: Produce bedgraphs for marking genes by proceeding to the subdirectory `bedgraph_visualizations/`
```
$ cd bedgraph_visualization
```

### STEP 4: Additionally, we would like to create bedgraphs by age and disease status, along with cell type.

This takes about 5 hours.
```
$ nohup python3 04_export_bigwig_by_cell_type_and_disease.py
```
