## Identify the GWAS SNPs overlapping cell-type specific peaks

### STEP 1: Iterate through all cell types

Run intersection of E2G linkages per cell type against the GWAS SNPs and store these results. Produce plot examining the different E2G linkages.

Run these scripts `01A_produce_enhancer_map.ipynb` and `01B_examine_E2G_links_per_cell_type.ipynb`.

### STEP 2: Visualize the E2G peaks that overlap with GWAS SNPs (or SNPs within close LD to these SNPs)

For each cell type, identify scE2G peaks that overlap with GWAS loci using `02A_run_intersection_analysis.py`. 

Then, visualize these results with `02B_R_examine_scE2G_overlap_results.ipynb`.