## Examine overlap between possible enhancer peaks in snATAC-seq and Spurrell enhancers

### STEP 0: Download Spurrell enhancers:

- There are from the Table S2A from Spurrell et al. 2022 has bulk ChIP-seq enhancer peaks for H3K27ac; this was saved as csv which we will convert to tsv bed format. There are 33317 enhancers called from the healthy donors.

### STEP 1: Extract the snATAC-seq peaks that are enhancers; make the Spurrell enhancers in bed format
# These enhancers were identified to be >1kb away from a TSS, so we will use this definition and the ChIPSeeker information using the interactive script
```
- 01_extract_enhancer_peaks.ipynb 
```
- This creates `01_snATAC_enhancer_peaks.bed`

### STEP 2: Perform degree of intersection between Spurrell enhancers and actual snATAC-seq enhancers, and count number of overlaps. 

Do this because the Spurrell enhancers are often quite long, with median length of 4244.0
```
$ bedtools intersect  -a Spurrell_healthy_enhancers.bed -b 01_snATAC_enhancer_peaks.bed -c >  enhancer_snATAC_peak_intersection.bed
```

Using `-c` counts the number of 500 nt peaks that overlap with the Spurrell enhancers. We are specifically interested in the number of Spurrell enhancers for which there are no intersecting snATAC-seq peaks, so we will use the $NF which includes the count of snATAC-seq peaks overlapping each of the 32525 
```
$ awk '$NF !=0 {print}' enhancer_snATAC_peak_intersection.bed | wc -l 
```

This gives 32525, which is number of 500 nt peaks that overlap with at least one of the Spurrell enhancers. This is 32525/33317 = 97.6% of the enhancers have at least 1 snATAC-seq peak.

### STEP 3: Perform null intersection to determine the expected number of Spurrell enhancers that will have at least 1 snATAC-seq peak overlapping it
```
$ nohup bash 03_examine_null_intersection.sh 
```

### STEP 4: Produce plots of the actual intersection vs the null distribution
```
- 04_examine_actual_enhancer_intersection_vs_null.ipynb
```

### STEP 5: Also examine the reciprocal situation, looking for the number of snATAC-seq peaks that have overlap with at least one of the Spurrell enhancers. 

In theory, since the snATAC-seq peaks are quite short, there shouldn't be many that overlap >1 Spurrell enhancer
```
$ bedtools intersect  -a 01_snATAC_enhancer_peaks.bed -b Spurrell_healthy_enhancers.bed -c >  snATAC_peak_to_enhancer_intersection.bed
$ awk '$NF !=0 {print}' snATAC_peak_to_enhancer_intersection.bed | wc -l 
```

The result is that 130291/566203 = 23% of these peaks overlap with Spurrell enhancers

# Using the ChIPSeeker annotations, we will see where these enhancers distribute
```
- 05_examine_peaks_overlapping_Spurrell_enhancers
```

It turns out that there is not a strong preference for any specific ChIPSeeker annotation (intron, exon, 3'UTR, 5'UTR, distal intergenic, downstream) for enhancer status.
