## Processing internal datasets

### STEP 0: Download the raw data from GEO and convert to fragment files

#### STEP 0A: Download the raw data from GEO

`00_directories.txt` is a list of the directories for the internal datasets, which are also listed below. Download from GEO the raw FASTQ files and move them to their respective sample directories within the `raw_data` directory. There should be 26 FASTQ subdirectories (one for each GEO accession), which correspond to 11 donors.

```
- Fetal-heart-8wk
- Fetal-LRV-1
- Fetal-LRV-14wk
- Fetal-LRV-2
- Fetal-LV-1
- Fetal-LV-2
- Fetal-OFT-1
- Fetal-OFT-14wk
- Fetal-OFT-2
- Fetal-wholeHeart-10wk
- H7CM-D20-1
- H7CM-D20-2
- H7CM-D328-1
- H7CM-D328-2
- K1727-LV-1
- K1727-LV-2
- LV-K1485-run-1
- LV-K1485-run-2
- LV-K1488-run-1
- LV-K1488-run-2
- LV-K1584-run-1
- LV-K1584-run-2
- LV-K1647-run-1
- LV-K1647-run-2
- LV-NDRI15755_10k
- LV-NDRI15755_5k 
```

#### STEP 0B: Run `cellranger-atac` for each directory
For each directory with fastq files, run the `cellranger-atac` pipeline. The version we used for `cellranger-atac` was 2.1.0.
```
nohup bash 00_send_cellranger_atac.sh &
```

### STEP 1: Copy over the fragment files in a new directory
The previous step will produce in each subdirectory a fragment file in count/outs/, which we will copy to a new directory called `fragment_files` for easier accessing.
```
$ nohup bash 01_get_fragment_files.sh & 
``` 

### STEP 2: Integrate with the donor metadata.
The metadata is `02_full_metadata.txt`. We will combine the metadata files with the external datasets together for the aggregated analysis (`../aggregated_analysis/`). The internal+external dataset metadata is located in `../aggregated_analysis/00_original_metadata.csv`.
