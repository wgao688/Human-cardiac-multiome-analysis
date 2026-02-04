## Identify the GWAS SNPs overlapping cell-type specific peaks

### STEP 1: Iterate through all cell types

Run intersection of E2G linkages per cell type against the GWAS SNPs and store these results. Produce plot examining the different E2G linkages.

Run these scripts `01A_produce_enhancer_map.ipynb` and `01B_examine_E2G_links_per_cell_type.ipynb`.

### STEP 2: Visualize the E2G peaks that overlap with GWAS SNPs (or SNPs within close LD to these SNPs)

For each cell type, identify scE2G peaks that overlap with GWAS loci using `02A_run_intersection_analysis.py`. 

Then, visualize these results with `02B_R_examine_scE2G_overlap_results.ipynb`.

### STEP 3: Intersect GWAS linkages with fetal reactivation

A: Identify fetal reactivation genes that have linkages in both our study (scE2G-based) and the original study using `03A_intersect_GWAS_with_fetalization_genes.ipynb`.
B: Visualize cell-type specific pseudobulked expression for DTL and ADAMTS7 using `03B_visualize_CM_RNA_fetalization_GWAS.ipynb`. This will also produce ggplot RNA pseudobulked plots.
C: Visualize cell-type specific pseudobulked accessibility for the peaks containing the linked SNPs using `03C_visualize_CM_ATAC_fetalization.ipynb`.
D: Produce ggplot ATAC pseudobulked plots using `03D_visualize_ATAC_violin.ipynb`.
