# Extracts relevant information in the R objects downloaded from GEO for Koenig et al. 2022 study so that they can be loaded in the adata format
library(tibble)
library(tidyr)
library(ggplot2)
library(Seurat)
library(Matrix)

load("GSE183852_DCM_Integrated.Robj")

# use ls() to identify which objects have been loaded
ls()

# create directory to store these inputs
adata_directory = "Seurat_to_h5ad/"
print(paste0("Creating directory called ", adata_directory))
dir.create(adata_directory)

# obtain the count matrix
# see here: https://www.youtube.com/watch?v=-MATf22tcak
print("Creating count matrix...")
count_matrix <- GetAssayData(RefMerge, assay = "RNA", slot = "counts")
writeMM(count_matrix, paste0(adata_directory, "matrix.mtx"))

# write the gene names
write.table(data.frame('gene'= rownames(count_matrix)),
            file = paste0(adata_directory, "gene_names.csv"), quote=F, row.names=F, col.names=F)

# metadata with the sex, names = cell type annotation, and the donor id, and condition (Donor vs. DCM), tech = SN vs. SC
relevant_metadata = RefMerge@meta.data[c("Sex", "Names", "orig.ident", "condition", "tech")]
write.csv(relevant_metadata, file = paste0(adata_directory, "metadata.csv"), row.names = T)

print("Script complete!")
