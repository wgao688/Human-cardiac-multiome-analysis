## Kanemaru 2023 spatial analysis

## RAW DATA (optional)

We ended up just using the processed data since the tiff files didn't work with spaceranger workflow, so this step can be skipped.

### STEP 1 (optional): Get metadata and fastq raw read paths. 

However, to download the raw data, `01_donor_info.txt` contains the donor metadata data and the fastq paths for wget. This info was manually fetched from [Array Express](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12975).

### STEP 2: Download the FASTQ files
```
$ nohup bash 02_download_fastqs.sh & 
```

Then, move these to `samples/` directory. 

### STEP 3: Manually fetch the links to the high resolution histology images.

They are located at [Heart Cell Atlas](https://www.heartcellatlas.org/). Copy the links for "High res image" into `03_image_metadata.txt`.

### STEP 4: Download the hi-res images based on the links in `03_image_metadata.txt`.
```
$ nohup bash 04_download_hires_images.sh
```

After this, I didn't to use the `.tiff` files in the `spaceranger` workflow, but the fiducials wouldn't align. Since we don't need the raw data, we can just proceed to analyzing the already processed data.

## PROCESSED DATA 

[Heart Cell Atlas](https://www.heartcellatlas.org/) already provides the raw counts and log normalized counts for the apical and LV samples. We will download those using their links, located in `processed/00_links.txt` using wget.

This produces four *.h5ad* files:
- `visium-OCT_AX_lognormalised.h5ad`: apex, log-normalized
- `visium-OCT_AX_raw.h5ad`: apex, raw
- `visium-OCT_LV_lognormalised.h5ad`: LV, log normalized
- `visium-OCT_LV_raw.h5ad`: LV, raw

For the fetal reactivation analysis, we focused on the LV samples. The analysis is located in `analysis_scripts/`.
