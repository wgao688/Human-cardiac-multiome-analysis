# Run the following to process the data

1. First,  `../raw_data` should contain all of the gzipped FASTQ files. Run this command to create separate directories for each of them
```
$ nohup bash 00_create_directories.sh -d ../raw_data_Dropseq/ & 
```

2. Then, run `STARSolo` on these fastq directory in a parallel fashion
```
$ nohup bash 01_run_STARSolo_parallel.sh -d ../raw_data_Dropseq/ &
```

3. Index the base files
```
$ nohup bash 02_index_bam_files.sh -d ../raw_data_Dropseq/ &
```

4. gzip compress the matrix files, as required by SoupX
```
$ nohup bash 03_gzip_STAR_files.sh -d ../raw_data_Dropseq/ &
```

5. Run SoupX using the count matrices in the raw/ and filtered/ for each sample subdirectory
```
$ nohup bash 04_send_SoupX.sh -d ../raw_data_Dropseq/ -s GeneFull &
```

```
6. Run SNP demultiplexing to confirm the individual donors
$ nohup bash 05_send_cellsnplite.sh -d ../raw_data_Dropseq/ &
```
