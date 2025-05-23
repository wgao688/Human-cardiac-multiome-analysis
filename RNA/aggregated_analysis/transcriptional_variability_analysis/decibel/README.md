## Perform transcriptional variability analysis using scallop. 

This will compute several different metrics for transcriptional noise, including Euclidean distance, Manhattan distance, and scallop noise. For further reading about the different noise metrics, please go to the [scallop paper](https://doi.org/10.7554/eLife.80380).

### STEP 0: Download packages for analysis; see LICENSE for this package

First, we need to download the decibel + scallop packages from the Ibanez-Sole et al. 2022.
```
$ conda create -n decibel
$ conda activate decibel
$ conda install python==3.7
$ pip install scallop
$ pip install jupyter
$ conda update setuptools
```

### STEP 1: Run decibel

Do this with subsampling of the counts and number of cells per cell type, since otherwise the decibel computation will take quite a while to perform. This will subsample the counts to 1000 UMIs per cell (and drop any cells for which there are fewer than 1000 UMIs). This is because the Euclidean distance and Manhattan distance metrics are sensitive to the number of UMIs. This will then perform the scallop analysis only for 100 cells per donor + cell type combination. 

Perform the subsampled version will take about 7-9 hours.

```
$ conda activate decibel
$ nohup python3 01A_run_decibel_with_subsampling.py &
```

#### If you wish to perform it on the cells (which will take about 60 hours), then perform the non-subsampled version.
```
$ conda activate decibel
$ nohup python3 01B_run_decibel_without_subsampling.py &
```

### STEP 2: Visualize the output of scallop interactively and produce final plots
```
- 02_analyze_decibel_subsampled_results.ipynb
```