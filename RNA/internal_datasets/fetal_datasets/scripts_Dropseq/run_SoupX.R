# Load the libraries
library(DropletUtils)
library(Seurat)
library(SoupX)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(Matrix)
library(ggrepel)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 1) {
  stop("Usage: Rscript run_SoupX.R  <input_directory: expected to be either Gene/ or GeneFull/>")
}

count_matrix_path <- args[1]

# load the raw and filtered matrices from STARSolo
print("Loading in the raw and filtered matrices. They should be gzipped already, otherwise this command will fail...")
raw.matrix <- Read10X(paste0(count_matrix_path, "/raw/"))
filt.matrix <- Read10X(paste0(count_matrix_path, "/filtered/"))

# create Seurat object of the filtered matrix
srat  <- CreateSeuratObject(counts = filt.matrix)
# create soup channel
soup.channel  <- SoupChannel(raw.matrix, filt.matrix)

# perform clustering
srat    <- SCTransform(srat, verbose = F)
srat    <- RunPCA(srat, verbose = F)
srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat    <- FindClusters(srat, verbose = T)

# add metadata
meta    <- srat@meta.data
umap    <- srat@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
soup.channel  <- autoEstCont(soup.channel)

# get the top 100 genes in ambient RNA
top_100_ambient_genes <- head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 100)
write.csv(top_100_ambient_genes, paste0(count_matrix_path, "top_100_ambient_genes.txt"))

# adjust the counts
adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

# get the mean number of counts per barcode removed
mean_counts_removed_per_cell = mean(colSums(adj.matrix) - colSums(filt.matrix))
write(mean_counts_removed_per_cell, file = paste0(count_matrix_path, "mean_ambient_RNA_removed.txt"))
saveRDS(adj.matrix, file = paste0(count_matrix_path, "soupX_corrected_counts.rds"))
