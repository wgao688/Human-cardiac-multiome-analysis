## Steps for preprocessing the Hill 2022 dataset

### STEP 1: Download based on the paper's data availability section

This study was focused on pediatric heart failure. Per the paper: "Raw and processed next-generation sequencing data have been deposited at the NCBI Gene Expression Omnibus under accession number GSE203275. The snRNA-seq data are available online at the Broad Single Cell Portal under study number SCP1852."

Therefore, to download these: log in to the portal (requires account). Then download `AllNuclei_snRNA_counts.mtx.gz`, `AllNuclei_snRNA_counts_colnames.txt.gz`, `AllNuclei_snRNA_counts_rownames.txt.gz`, `AllNuclei_snRNA_metadata.csv`

Then gunzip these files
```
gunzip *.gz
```

### STEP 2: Perform iteractive processing using the jupyter notebooks
- `01_load_Hill_to_adata.ipynb`: This will load in the metadata and count files to create anndata
- `02A_preprocess_Hill_postnatal_ND_adata.ipynb`: Process the non-diseased donors
- `02B_preprocess_Hill_postnatal_diseased_adata.ipynb`: Process the diseased donors

### OUTPUTS:
- `02_processed_Hill_ND.h5ad`: Processed non-diseased LV adata
- `02_processed_Hill_D.h5ad`: Processed diseased LV adata
