### This is the directory to perform trajectory analyses from ND to diseased cells 
### Perform this for Cardiomyocyte, Fibroblast, Endothelial, Myeloid

# STEP 1: send the trajectory script. 
Run `bash 01_send_trajectory.sh`, which will iteratively run `01_run_trajectory_per_cell_type.py` for all cell types. This python script uses functions that are in `trajectory_funcs.py`. The trajectory is calculated according to the density of the cells in each state in the embedding space.

This takes about 1 hour for each cell type.

The results are the visualization of the trajectories from fetal->ND and from ND->DCM/HCM/ICM for each of the major cell types (cardiomyocyte, endothelial, fibroblast, myeloid). All of the genes are hierarchically clustered along the density gradient. 

### STEP 2: Perform ORA for hierarchical clusters

For the hierarchical clusters identified for each trajectory, we will perform ORA to identify pathways enriched.