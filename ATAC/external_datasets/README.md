## Downloading fragment files

This directory stores all of the fragment files for the each of the datasets, and their indices
 
- `ENCODE`: The ENCODE fragments files were downloaded from the ENCODE v4 portal
- `Kanemaru`: The Kanemaru fragment files were downloaded from the online portal, under ATAC fragments: https://www.heartcellatlas.org/. Only keep the LV samples.
- `Kuppe`: The Kuppe fragment files were downloaded from the Zenodo archive listed in the "Data Availability" section: https://zenodo.org/records/6578047
- `Ameen`: The Ameen datasets were downloaded from GEO accession GSE181346
- Penn: The Penn fragment files were generated in `../internal_datasets` using the raw FASTQ files. Download these from this study's GEO.

All further preprocessing (count matrix generation, QC analysis, intergration, etc.)  will be performed together in `../aggregated_analysis/`
