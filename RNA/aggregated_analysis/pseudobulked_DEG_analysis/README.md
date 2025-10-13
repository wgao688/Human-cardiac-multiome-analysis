## Differential gene expression analysis

This directory contains the scripts for the DEG analysis across sex, aging, development and disease.

### STEP 1: Determine donors with enough UMIs for stable estimation

Since the number of cells per cell type and donor varies significantly, we will determine how many UMIs are needed to obtain a representative transcriptome profiling per donor. We will do this by subsampling the number of reads for each cell-type specific donor pseudobulked profile, using the 5 donors with the most number of reads per cell type. Using the script  `01_determine_number_cells_for_pseudobulk_analysis.ipynb`, we estimated that >50K UMIs in the pseudobulked profile will provide a spearman correlation >0.8 to the full profile.

### STEP 2: Generate the pseudobulked counts.

Now, we will use `02_create_pseudobulk_raw_counts.ipynb` to generate the donor pseudobulked cell-type specific count matrices and metadata for downstream DEG analysis.


### Subdirectories

The subdirectories include the following:
- `fully_pseudobulked_analysis/`: this includes the analysis where we compared the "batch as cov" approach and the "meta-analysis" against the silver standard bulk RNA-seq left ventricle GTEx data.
- `saturation_analysis`: this includes analysis where we downsample the number of donors to see how the number of DEGs changes as a function of donors, when using the "batch as cov" approach.
- `sex_and_aging_DEG_analysis`: scripts for identifying sex (male vs. female) and aging (old vs young) DEGs using "batch as cov" and "meta-analysis" approaches
- `developmental_DEG_analysis`: scripts for identifying developmental DEGs (fetal vs. young)
- `disease_DEG_analysis`: scripts for identifying disease DEGs, either binarized (Y vs N) or non-binarized (DCM, HCM, and ICM vs. non-diseased)
- `intersection_analysis`: this includes scripts for identifying the intersection of DEGs across contrasts. This includes scripts for further examination fetal reactivation in diseased state.
