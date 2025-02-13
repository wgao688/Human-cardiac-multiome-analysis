# Perform cell type proportion analysis using scanpro

# STEP 1: Cell type proportion analysis. 

Process the adata and produce plots for raw cell type proportions. We will also fit a model that takes into account the relevant covariates (age, disease status), but does not perform extra steps related to the interdependency and mean-variance relationship which propeller and scanpro do. This script will also produce data files for scanpro. Perform this interactively.
```
- 01_setup_cell_type_proportion_analysis.ipynb
```

# STEP 2: Run scanpro analysis interactively

This takes into account the covariates and the mean-variance relationship (runs propeller in python)
```
- 02_run_scanpro.ipynb 
```

# STEP 3: Perform plots using ggplot
```
- 03_produce_cell_type_proportion_plots.ipynb
```
