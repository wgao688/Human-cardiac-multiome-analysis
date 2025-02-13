### ArchR analysis

We will load in the processed ATAC-seq dataset into ArchR for applications that require usage of ArchR, such as the useful genome-browser bigwig options and peak-gene correlation.

### STEP 1: Convert from SnapATAC2 format to ArchR compatible format

#### STEP 1A: Export the fragments to bed files. This takes about 1.5 hours
```
$ nohup python3 01B_export_fragments.py &
```

#### STEP 1B: Tabix the fragment files
```
$ nohup bash 01B_tabix.sh & 
```

### STEP 2: Run ArchR in this directory: `run_ArchR/`
