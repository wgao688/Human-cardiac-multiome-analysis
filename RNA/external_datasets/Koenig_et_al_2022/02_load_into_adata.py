# Step 01 was an R script that converted the data from the Seurat object into a format can be read into adata format
# Step 02 (this script) will generate the adata file

import pandas as pd
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import matplotlib.pyplot as plt

# read in the sparse matrix, and add it to a new adata object
X = io.mmread("Seurat_to_h5ad/matrix.mtx")
adata = anndata.AnnData(X = X.transpose().tocsr())

# load in the metadata
# index_col=0 reads the first column (rownames in R) as index
metadata = pd.read_csv("Seurat_to_h5ad/metadata.csv", index_col = 0)

# load in the genes line by line
with open("Seurat_to_h5ad/gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

adata.obs = metadata
adata.var.index = gene_names

adata.write("Koenig_all_sc_snRNA.h5ad")
print("Script complete!")
