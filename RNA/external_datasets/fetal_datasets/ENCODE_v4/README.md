## Steps for preprocessing the ENCODE v4 dataset (fetal donors)

### STEP 0: Download data
I manually generated the description of each of the datasets in `00_fetal_multiome_metadata.txt` using the link provided in `https://www.encodeproject.org/search/?type=MultiomicsSeries&lab.title=Michael+Snyder%2C+Stanford&related_datasets.replicates.library.biosample.life_stage=embryonic&biosample_ontology.term_name=heart` as there are fewer than 10 donors.

The following ENCODE accessions were analyzed:
- ENCFF069ATM 
- ENCFF248EWR
- ENCFF684YRB
- ENCFF727JRO
- ENCFF775ANN
- ENCFF776DQR
- ENCFF802AQC
- ENCFF805YRY
- ENCFF849ALE

Then, use the `00_download_tarballs.sh` to download the count matrix tarballs using wget, which fetches the fastq files. Each of those ENCODE accessions will have its own directory with the relevant files after running `00_download_tarballs.sh`.

### STEP 1: Reformat the data to extract the count matrices
```
# unpack the tarballs
$ bash 01A_unpack_tarballs.sh

# gzip the extracted files
$ bash 01B_gunzip_files.sh
```

### STEP 2: Generate the combined and processed adata from these unfiltered matrices interactively
- `02_generate_combined_adata.ipynb`

This will produce an intermediate adata file called `02_fetal_ENCODE.h5ad`.

### STEP 3: For the combined adata, perform cell type annotation
Similar to the postnatal datasets, the fetal datasets have not been given cell type annotations. Use the interactive jupyter notebook `03_annotate_adata.ipynb` to annotate the cell types. 

We will then reformat the adata using `03_reformat_adata.ipynb`. This will produce the final adata called `03_reformatted_ENCODE_fetal.h5ad`. 

### OUTPUTS:
- `03_reformatted_ENCODE_fetal.h5ad`: This is the combined adata containing whole heart fetal samples.
