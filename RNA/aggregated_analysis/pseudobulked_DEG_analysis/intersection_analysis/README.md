### Perform the intersection analysis between DEG sets

### STEP 1: Perform intersection analysis
Check to see if there DEG overlap in concordant directions is statistically significant or not, given the set of genes. Perform using this script `01_intersect_analysis.py`, which calls functions from `intersection_scripts.py`

### STEP 2: Visualize results

Using this script `02_produce_intersection_plot.ipynb`

### STEP 3: Perform over-representation analysis for fetal reactivation genes

Identify the genes overlapping in concordant up or down for disease and fetal, and perform over-representation analysis (ORA) using the Hallmark pathways using this script `03A_perform_ORA_for_fetalization.ipynb`. Visualize results using this script `03B_visualize_top_fetalization_ORA_R.ipynb`.
