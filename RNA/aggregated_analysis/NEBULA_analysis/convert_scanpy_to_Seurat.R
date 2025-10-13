# Converts scanpy anndata to Seurat object
# Adapted from helpful video: https://www.youtube.com/watch?v=4zqMV1N96sY

# Load necessary libraries
library(reticulate)
library(Seurat)
library(argparse)
library(Matrix)
start_time <- Sys.time()

# Create a parser
parser <- ArgumentParser(description = 'Convert scanpy annData to seurat object')

# Add arguments
parser$add_argument('-a', '--adata_path', required=TRUE, help='Path to the input .h5ad file')
parser$add_argument('-o', '--output_directory', required=TRUE, help='Directory to save the output files')

# Parse the arguments
args <- parser$parse_args()

print("Reading in the .h5ad file in Python and converting it to the directory format for Read10X() in Seurat")

# Assign arguments to variables
adata_path <- args$adata_path
output_directory <- paste0(args$output_directory, "/")

# Print the paths to confirm
cat("Using the following paths:\n")
cat("adata_path: ", adata_path, "\n")
cat("output_directory: ", output_directory, "\n")

# Use Python in R with reticulate
use_python("/home/william/anaconda3/envs/r-kernel/bin/python3")  # Modify the path to your Python executable if necessary

# Create the output directory if it doesn't exist
dir.create(output_directory, showWarnings = FALSE)

# Python code executed within R
# Python code executed within R
py_run_string(sprintf("
import scanpy as sc
from scipy import io
import os
import gzip
import shutil
import pandas as pd

adata_path = '%s'
output_directory = '%s'

# Read H5AD file using scanpy
adata = sc.read_h5ad(adata_path)

# use the raw counts, make sure they are integers
adata.X = adata.layers['counts'].astype(int)

print(adata)

os.makedirs(output_directory, exist_ok=True)

# Create barcode file in output directory
with open(output_directory + 'barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\\n')

# Create the features file
with open(output_directory + 'features.tsv', 'w') as f:
    for item in ['\\t'.join([x, x, 'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\\n')

# Create the matrix file - this can take quite a while
io.mmwrite(output_directory + 'matrix.mtx', adata.X.T)

# Gzip the files
with open(output_directory + 'barcodes.tsv', 'rb') as f_in:
    with gzip.open(output_directory + 'barcodes.tsv.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

with open(output_directory + 'features.tsv', 'rb') as f_in:
    with gzip.open(output_directory + 'features.tsv.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

with open(output_directory + 'matrix.mtx', 'rb') as f_in:
    with gzip.open(output_directory + 'matrix.mtx.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

# Remove the original uncompressed files
os.remove(output_directory + 'barcodes.tsv')
os.remove(output_directory + 'features.tsv')
os.remove(output_directory + 'matrix.mtx')

# Save UMAP embeddings if they exist
if 'X_umap' in adata.obsm.keys():
    pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names).to_csv(output_directory + 'umap_embeddings.csv')

# Transfer the metadata annotations over as well
adata.obs.to_csv(output_directory + 'metadata.csv')
", adata_path, output_directory))

print("Now creating the Seurat object in R")

# Continue in R: Read the matrix and metadata in Seurat
object_counts <- Read10X(output_directory)
#object_counts <- readMM(ouput_directory)

# Read in the metadata
scanpy_metadata <- read.csv(file.path(output_directory, "metadata.csv"))
rownames(scanpy_metadata) <- scanpy_metadata[["X"]]
scanpy_metadata[["X"]] <- NULL  # remove this column

# Now create the Seurat object with the raw count matrix and metadata
obj <- CreateSeuratObject(counts = object_counts, meta.data = scanpy_metadata)

# If UMAP embeddings were saved, add them to the Seurat object
umap_file <- file.path(output_directory, "umap_embeddings.csv")
if (file.exists(umap_file)) {
  umap_embeddings <- read.csv(umap_file)
  umap_embedding_matrix <- as(umap_embeddings[, 2:3], "matrix")
  rownames(umap_embedding_matrix) <- umap_embeddings[, 1]
  colnames(umap_embedding_matrix) <- c("UMAP1", "UMAP2")
  obj[["umap"]] <- CreateDimReducObject(embeddings = umap_embedding_matrix, key = "UMAP_", global = T, assay = "RNA")
}

# Save the Seurat object
saveRDS(obj, file = paste0(output_directory, "seurat_object.rds"))
end_time = Sys.time()
execution_time <- end_time - start_time
cat("Execution time:", execution_time, "\n")
print("Script complete!")
