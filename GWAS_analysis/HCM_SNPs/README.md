## HCM GWAS SNP analysis

### STEP 0: Download the GWAS hits

We will examine the GWAS SNPs from Tadros et al. 2025 [Large-scale genome-wide association analyses identify novel genetic loci and mechanisms in hypertrophic cardiomyopathy](https://www.nature.com/articles/s41588-025-02087-4).

The SNPs are from Supplemental Table 8.

### STEP 1A: Convert the Tadros SNPs, which are in GRCh37 (hg19) coordinates to hg38 coordinates; downloaded from UCSC the hg19 to hg38 liftOver chain. 

Run `liftOver`:
Use `01A_Tadros_visualize_GWAS_metadata.ipynb` to inspect the metadata for the Tadros 68 SNPs and save to a bed-like format. 
But to run this with `liftOver`, we need to delete the column titles. 

```
$ liftOver Tadros_GWAS_hits_GRCh37.bed hg19ToHg38.over.chain.gz Tadros_GWAS_hits_GRCh38.bed Tadros_GWAS_hits_GRCh37.unmapped.bed
```

- `Tadros_GWAS_hits_GRCh38.bed` now contains the coordinates for these SNPs in hg38

### STEP 2: Create list of SNPs for PLINK, and also create a bed file format of the combined SNPs

Use this script: `02_create_SNP_list.ipynb`
