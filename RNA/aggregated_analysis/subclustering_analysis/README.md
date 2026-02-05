## Subclustering analysis

In this subdirectory, we will perform subclustering at the cell type level for the final RNA adata. We will try to see if there are specific subclusters of each cell type shared across studies. The subclusters will be analyzed using cNMF and Silhouette analysis to determine the best leiden cluster resolution.

### STEP 1: Run subclustering to each cell type, using 8 different leiden resolutions: 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0.

A: Create the adata files that are subsetted to each cell type. These will be stored here: `pre_subclustering/`
```
$ nohup python3 01A_create_cell_type_specific_adata.py &
```

B: This will send `01B_run_subcluster_for_one_cell_type.py` non-interactively for each cell type. That script will subset to the specified cell-type and then perform leiden clustering at the 8 resolutions specified. Since this waits 30 minutes between each cell type (since scVI is computationally expensive, this will take around 8 hours to run). 
```
nohup bash 01B_send_subclustering.sh & 
```

### STEP 3: Perform silhouette analysis

This script will iterate through each cell type, to identify the best resolution based on the maximum silhouette score.
```
nohup python3 03_run_silhouette_single_cell_type.py -c "cell_type" & 
```

### STEP 3: Perform manual subclustering

Now, we will use the appropriate leiden cluster resolution and marker genes / literature search to perform manual subclustering. The scripts for this are in the subdirectory `perform_annotation/`.
