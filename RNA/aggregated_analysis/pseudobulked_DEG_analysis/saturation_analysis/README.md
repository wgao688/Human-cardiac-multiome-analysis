## DEG saturation analysis

Examine how increasing the number of donors affects the number of DEGs identified, by subsampling the number of donors for cardiomyocytes and fibroblasts.


### STEP 1
Perform the subsampling and DE analysis
```
nohup python3 01A_run_Cardiomyocyte_saturation_analysis.py & 
nohup python3 01B_run_Cardiomyocyte_saturation_analysis.py &
```

### STEP 2
Visualize in R
```
- 02_plot_saturation.ipynb
```
