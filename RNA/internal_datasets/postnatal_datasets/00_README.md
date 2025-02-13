### This directory performs all preprocessing for postnatal Penn donors

The directories here include `scripts/` for the non-interactive scripts, `combined_analysis/` for the jupyter notebooks containing the scripts for downstream analysis after combining all of the count matrices together. 

### STEP 1: Download the data

Please download the FASTQ files from GEO. `01_metadata.txt` specifies the donor information for each of these libraries. There should be 54 pairs of FASTQ R1 and R2 downloaded, including all of these:

- K1430-LV
- K1430-LV-2nd
- K1430-LV-FACS
- K1430-LV-FACS-KAPA
- K1485-LV-FACS
- K1485-LV-nonFACS
- K1488-LV-FACS
- K1488-LV-nonFACS
- K1545-LV-FACS
- K1545-LV-run123
- K1549-LV-FACS
- K1549-LV-nonFACS
- K1561-LV-FACS
- K1570-LV-FACS
- K1570-LV-nonFACS
- K1584-LV-FACS
- K1584-LV-FACS-e2
- K1584-LV-run1n2
- K1617-LV
- K1617-LV-2nd
- K1617-LV-FACS
- K1617-LV-FACS-KAPA
- K1622-LV-FACS
- K1622-LV-run123
- K1647-LV-FACS
- K1647-LV-nonFACS
- K1702-LV-RE-FACS
- K1702-LV-RE-nonFACS
- K1718-LV-FACS
- K1718-LV-nonFACS
- K1727-LV-FACS
- K1727-LV-nonFACS
- K1739-LV-FACS
- K1739-LV-nonFACS
- K1783-LV-FACS
- K1783-LV-nonFACS
- K1783-LV-RE-FACS
- K1783-LV-RE-nonFACS
- ND15755-Atr-1st
- ND15755-Atr-3rd
- ND15755-LV-1st-run123
- ND15755-LV-2nd-run1234
- ND15755-RV-1st
- ND15755-RV-2nd
- ND15755-RV-3rd
- ND15840-LV-1st-run1n2
- ND15840-LV-2nd-run1n2
- ND15840-LV-3rd
- ND15840-RA-2nd-run1n2
- ND15840-RA-3rd
- ND15840-RA-4th
- ND15840-RV-1st
- ND15840-RV-2nd
- ND15840-RV-3rd

### STEP 2: Create raw data directory. 

Place all of the postnatal FASTQ files in the same directory called `raw_data/`. Then, run `00_create_directories.sh`, which will create one directory for each fastq pair, based on the prefix in the FASTQ file names:

```
$ nohup bash 00_create_directories.sh -d raw_data &
```

### STEP 3: Run STARSolo to align to the reference genome
With these fastq files in their respective directories, `01_run_STARSolo_parallel.sh` will run STARSolo in a parallel fashion on the Dropseq data, resulting in BAM files and count matrices. Move these into the `raw_data/` directory.

### STEP 4: Index the bam files and gzip the count matrix directories:
```
# index the bam files
02_index_bam_files.sh

# compress the count matrix directories
03_gzip_STAR_files.sh
```

### STEP 5: Run SoupX in a parallel fashion (5 at a time): 
```
nohup bash scripts/send_SoupX.sh -d raw_data/ -s GeneFull -p 5 &
```

### STEP 6: Combine the SoupX corrected count matrices in .RDS file into a format that can be loaded as an adata in scanpy

First, combine all of the SoupX count matrices together using this script
- `combined_analysis/01_combine_SoupX_corrected_count_matrices.ipynb`

Then, convert this from an RDS object to .h5ad object
- `combined_analysis/02_scanpy_combine_adata.ipynb`

Run the standard preprocessing that we use for both internal and external datasets, interactively:
- `combined_analysis/03_preprocess_Penn_diseased_LV.ipynb`

Lastly, perform separate filtering for diseased and ND donors:
- `combined_analysis/03_preprocess_Penn_diseased_LV.ipynb`
- `combined_analysis/03_preprocess_Penn_non_diseased_LV.ipynb`
