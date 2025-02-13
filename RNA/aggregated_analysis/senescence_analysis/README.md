## Senescence analysis

In this directory, we will explore whether diseased or aged cardiac cells displayed increased senescence. We will do this using the `scanpy.tl.score_genes()` function in scanpy that inputs an adata object and a list of genes. This gene score then represents the average expression of the gene set compared to a "refernence" set of genes belonging to the same expression bin.

### STEP 0: Download gene sets:
- The reactome senescence gene set: R-HSA-2559583, which the SenMayo paper compared their gene set against. Use the script to download this gene set and format it for scanpy: `00_extract_Senescence_gene_set_list`
- The SenMayo gene set is here: `SenMayo_list.txt`. I extracted these genes from the Saul et al. 2022 supplemental materials.
- A third list will just be the two senescence markers CDKN1A (p16) and CDKN2A (p21).

### STEP 1: Compute the gene scores at single cell type. 

We will use the scanpy function to get the scores for each cell type, and then pseudobulk to obtain a mean senescence score at the cell type and donor-level, so that we can perform covariate modeling to account for batch effects. Perform this interactively with the jupyter notebook. 
``
- 01_senescence_analysis.ipynb
```

# STEP 2: Produce ggplots for the senescence analysis in R.
```
- 02_produce_senescence_analysis_plots.ipynb
```
