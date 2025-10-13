# Download the hg38 VCF files from the 1000 genomes project, identify SNPs in high LD (r^2 > 0.8) to the GWAS SNPs 

The 1000 genomes LD files for plink are already in the `GWAS_DCM` directory, so we will skip the first two steps which address how to obtain this file and interpret PLINK results.

### STEP 1: For each HCM SNP, query around 100kb for SNPs and compute r^2 using the 1000 genomes project. This takes about 10 minutes.
```
$ nohup bash 01_obtain_LD_results.sh &
```
--> all SNPs except one: 11:47343463 had a corresponding rsID in the 1000 genomes project

### STEP 2: For each SNP, identify other SNPs that have r^2>0.2
```
- 02_identify_high_LD_SNPs.ipynb
```

### STEP 3: Perform analysis for using the 1000 genomes project, European donors only

```
- 03_obtain_LD_EUR_results.sh
```

### STEP 4: Identify SNPs that have r^2>0.2
