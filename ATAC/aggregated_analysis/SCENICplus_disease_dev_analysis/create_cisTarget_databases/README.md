## Create custom cistarget db

This directory was largely copied over from the SCENIC+ tutorial, but need to edit the following files. Make sure to use the scenicplus conda environment, an environment that has SCENIC+ dependencies installed. See here: https://scenicplus.readthedocs.io/en/latest/install.html

This follows the tutorial for [Creating custom cistarget database](https://scenicplus.readthedocs.io/en/latest/human_cerebellum_ctx_db.html).

### STEP 1: Create a fasta with the consensus peak sequences for the cardiac dataset, using hg38.fa
```
- nohup bash 01_create_fasta_from_consensus.sh &
- The output of this is human_cardiac_peaks.fa
```

### STEP 2: Using these sequences, create a cistarget database. This takes a while (about 11 hours) to run.
```
- nohup bash 02_create_cistarget_motif_db.sh &
```

### OUTPUTs: custom cistarget databases for SCENIC+
```
- human_cardiac.motifs_vs_regions.scores.feather
- human_cardiac.regions_vs_motifs.rankings.feather
- human_cardiac.regions_vs_motifs.scores.feather
```
