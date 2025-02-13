### STEP 1: Download data

- Please download the raw FASTQ files from GEO, and then save them to the directory `raw_data_Dropseq`, which should have these directories, each containing the R1 and R2 FASTQ files:

- `Fetal-1st-LV-0315-1-run1n2`
- `Fetal-1st-LV-0315-2-run1n2`
- `Fetal-2nd-Atria-1`
- `Fetal-2nd-Atria-2`
- `Fetal-2nd-OFT-1`
- `Fetal-2nd-OFT-2`
- `Fetal-Atria-18wk`
- `Fetal-LRV-18wk-male1-run1n2`
- `Fetal-LRV-18wk-male2-run1n2`
- `Fetal-LRV-18wk-male3-run1n2`
- `Fetal-LRV-18wk-male4-run1n2`
- `Fetal_LV_18wk_e1-run1n2`
- `Fetal-LV-18wk-e2-run1n2`
- `Fetal-LV-18wk-e3`

### STEP 2: Perform analysis. 

Do the analysis in the directory called `combined_analysis` according to its `README.md` file. That directory will call scripts from `scripts_Dropseq`.

### STEP 3: Updated metadata.

To verify the donors for each library, we ran cellsnplite as part of the analysis in `combined_analysis`. Run `01_view_cellsnplite_results.ipynb` to visualize the SNP calling results. After this step, we will store the metadata for these fetal samples here in `01_fetal_updated_metadata.csv`
