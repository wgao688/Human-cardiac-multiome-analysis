# Human cardiac multiome analysis
This repository contains the scripts used for the analysis of the human cardiac snRNA-seq, snATAC-seq, and spatial transcriptomic datasets in the manuscript titled "Single-cell multiomic integration identifies widespread, cell-type resolved 
fetal reactivation in the diseased human heart"). For any questions regarding these scripts, please reach out to william.gao@pennmedicine.upenn.edu. 

This repository includes several subdirectories. For each subdirectory, there is a `00_info.txt` file describing the contents of the directory, and the order of the commands used to generate the paper figures. Some of the scripts are run in the command line; others are interactive Jupyter notebooks. To run these notebooks/scripts, there are some R/python libraries that need to be downloaded, which will be indicated in the scripts/notebooks. Because some files are too large, they are not included in this directory, but can be regenerated using the scripts, as described in the `00_info.txt` files within each directory. These large files include *.h5ad (anndata format) objects with snRNA-seq + snATAC-seq + spatial objects, and FASTQ files.

The new samples generated as part of this study (19 snRNA-seq donors, 11 snATAC-seq donors) are deposited as raw FASTQ files in the Gene Expression Omnibus under the accession number *TBD*. 
