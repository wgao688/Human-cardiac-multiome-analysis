## DCM GWAS SNP analysis

### STEP 0: Download the GWAS hits

We will examine the GWAS hits from Jurgens et al. 2024 [Genome-wide association study reveals mechanisms underlying dilated cardiomyopathy and myocardial resilience](https://www-nature-com.proxy.library.upenn.edu/articles/s41588-024-01975-5#Sec24) and Zheng et al. 2024 [Genome-wide association analysis provides insights into the molecular etiology of dilated cardiomyopathy](https://www-nature-com.proxy.library.upenn.edu/articles/s41588-024-01952-y#Sec32).

The SNPs are extracted the information of interest from Supplemental Table 40 of Jurgens et al. 2024, which performed a comparison between their study and Zheng et al. 2024.

### STEP 1A: Convert the Jurgens SNPs, which are in GRCh37 (hg19) coordinates to hg38 coordinates; downloaded from UCSC the hg19 to hg38 liftOver chain. Run liftOver:

Use `01A_Jurgens_visualize_GWAS_metadata.ipynb` to inspect the metadata for the Jurgens 65 SNPs and save to a bed-like format. 
But to run this with `liftOver`, we need to delete the column titles. 

```
$ liftOver Jurgens_GWAS_hits_GRCh37.bed hg19ToHg38.over.chain.gz Jurgens_GWAS_hits_GRCh38.bed Jurgens_GWAS_hits_GRCh37.unmapped.bed
```

- `Jurgens_GWAS_hits_GRCh38.bed` now contains the coordinates for these SNPs in hg38

### STEP 1B: Convert the Zheng SNPs, which are also in GRCh37 (hg19) coordinates to hg38 coordinates

Use `01B_Zheng_visualize_GWAS_metadata.ipynb` to obtain the bed-like file, and then run `liftOver`
```
$ liftOver Zheng_GWAS_hits_GRCh37.bed hg19ToHg38.over.chain.gz  Zheng_GWAS_hits_GRCh38.bed Zheng_GWAS_hits_GRCh37.unmapped.bed
```

### STEP 2: Examine the intersection between the two GWAS studies in terms of SNPs

Use this script: `02_examine_hg38_GWAS_intersection.ipynb`

### STEP 3: Create combined list of SNPs for PLINK, and also create a bed file format of the combined SNPs

Use this script: `03_create_SNP_list.ipynb`
