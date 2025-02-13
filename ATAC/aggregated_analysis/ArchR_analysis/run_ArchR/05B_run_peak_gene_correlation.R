# Add gene expression (from sce) to the ArchR project, saving a new ArchR object
# Then perform peak-gene correlation between all genes within 1Mb of the gene

library(ArchR)
library("unixtools")
set.tempdir("/mnt/data1/william/tmp")
addArchRThreads(threads = 32)
addArchRGenome("hg38")
library(pheatmap)
`%notin%` = Negate(`%in%`)

#library(Seurat)
library(SingleCellExperiment)

start_time = Sys.time()
working_dir <- getwd()

proj3 <- loadArchRProject(path = working_dir,
                          force = FALSE,
                          showLogo = FALSE)

metadata <- getCellColData(proj3)

# load in the sce with the granges
sce <- readRDS("../make_Multiome_sce/final_RNA_sce_with_granges.rds")

# get the ArchR indices found in the RNA modality
colnames(sce) <- sce$ArchR_index

# filter proj3 (ATAC modality) to just the multiome nuclei from in the RNA modality 
filt_proj3 <- proj3[ sce$ArchR_index, ]

# make the rowranges consistent with "chr" added
rr <- rowRanges(sce)
new_seqnames <- paste0("chr", seqnames(rr))
# update seqlevels to match new seqnames
seqlevels(rr) <- unique(new_seqnames)
# assign the modified seqnames back
seqnames(rr) <- new_seqnames
# update the rowRanges in the SCE object
rowRanges(sce) <- rr
# also add the genes as rownames
rownames(sce) <- rowRanges(sce)$gene_name

# add RNA expression to ArchR project
filt_proj3 <- addGeneExpressionMatrix(
  input = filt_proj3, 
  seRNA = sce, 
)

# add peak matrix
filt_proj3 <- addPeakMatrix(filt_proj3)

# perform peak-to-gene linkage
filt_proj3 <- addPeak2GeneLinks(
  ArchRProj = filt_proj3,           
  useMatrix = "GeneExpressionMatrix",  # matrix for gene-peak linkage
)

# save the ArchR project in a new directory
saveArchRProject(ArchRProj = filt_proj3)

end_time = Sys.time()

time_taken <- as.numeric(difftime(end_time, start_time, units = "hours"))
cat("Total time taken:", time_taken, "hours\n")
