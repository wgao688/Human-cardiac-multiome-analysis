## SCENIC+ analysis

In this directory, we will run SCENIC+ analysis.
Make sure that you have SCENIC+ installed in a separate conda environment, as it requires more specific dependencies. See here: https://scenicplus.readthedocs.io/en/latest/install.html for installation instructions.

### STEP 1: Generate subsamples of the ATAC and RNA adata by cell types and age+disease status. 

We will only perform analysis for cell types in which there are at least 750 cells for all age+disease statuses (fetal, ND, disease). These are the following cell types that satisfy those conditions: Cardiomyocyte, Endothelial, Fibroblast, Lymphoid, Myeloid, Pericyte
```
- 01A_generate_subsampled_ATAC.ipynb
- 01B_generate_subsampled_RNA.ipynb
```

### STEP 2: Prepare the snapatac2 adata for pycistopic
```
- 02_prepare_snapatac2_for_pycistopic.py
```

### STEP 3A: Run pycistopic. This takes about 5 hours.
```
$ conda activate scenicplus
$ nohup python3 03A_run_pycistopic.py &
```
### STEP 3B: Create the cistarget database. Follow steps in that directory.
```
$ cd create_cisTarget_databases/
```

### STEP 4: Run SCENIC+, copying over the scripts
Run this in an environment that has SCENIC+ downloaded according to its installation instructions. 
```
$ conda activate scenicplus
$ scenicplus init_snakemake --out_dir scplus_pipeline
$ cd scplus_pipeline/Snakemake/
```

Before running, we will need to custom the snakemake config file (config.yaml). Edit the paths accordingly.
```
$ snakemake --cores 40
```

After this step is done, we will extract the activity scores for TFs across the different states (fetal, ND, D) using these notebooks:
```
- 04A_analyze_scenicplus_outs.ipynb
- 04B_ggplot_SCENIC_regulons.ipynb
```

### STEP 5: Examine some fetalization TFs and the TF->enhancer->gene linkages
```
- 05_SCENICplus_triplet_overlap_with_fetalization.ipynb
```
