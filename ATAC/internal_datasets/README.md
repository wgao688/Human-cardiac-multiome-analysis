## Processing internal datasets

`00_directories.txt` is a list of the directories for the internal datasets. Download from GEO the raw FASTQ files and move them to their respective sample directories.

### STEP 0: Run cellranger-arc for each directory
For each directory with fastq files, run the cellranger-atac pipeline.
```
nohup bash 00_send_cellranger_atac.sh &
```

### STEP 1: Copy over the fragment files in a new directory
```
$ nohup bash 01_get_fragment_files.sh & 
```

### STEP 2: Examine the donor metadata. 
The metadata in `02_full_metadata.txt`.

### STEP 3: Run cellsnplite on the Fetal samples to confirm their identities.
Call the script `run_cellsnplite_pseudobulked.sh` for each of the FASTQ directories.
```
$ bash 03_send_cellsnplite.sh -d .
```

### STEP 4: Verify the sex of all of the donors
`- 04_verify_sex_BAM.ipynb`
