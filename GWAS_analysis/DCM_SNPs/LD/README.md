# Download the hg38 VCF files from the 1000 genomes project, identify SNPs in high LD (r^2 > 0.8) to the GWAS SNPs 

### STEP 1: Download the high coverage 1000 genomes project VCF files
```
$ nohup bash 01_download_vcf.sh &
```

### STEP 2: Read the plink notes for more detail on how to use PLINK, located in `02_plink_notes.txt`
# First, we will prepare the plink output files, one for each chromosome. This takes about 2 hours.
```
$ nohup 02_prepare_plink_files.sh & 
```

### STEP 3: Then, for each SNP, query around 1 megabase for SNPs and compute r^2 using the 1000 genomes project. This takes about 10 minutes.
```
$ nohup bash 03_obtain_LD_results.sh &
```
--> all SNPs except one: 11:47343463 had a corresponding rsID in the 1000 genomes project

### STEP 4: For each SNP, identify other SNPs that have r^2>0.8
```
- 04_identify_high_LD_SNPs.ipynb
```
### STEP 5: Extract LD blocks for European ancestry only

A: Extract the European donors in 1K genomes project
B: Extract the VCF for the European ancestry only

```
- nohup bash 05B_extract_EUR_vcf.sh & 
```
