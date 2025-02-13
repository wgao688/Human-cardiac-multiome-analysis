# Performs QC filtered, batch correction, and harmonization

library(ArchR)
library("unixtools")
set.tempdir("/mnt/data1/william/tmp")
addArchRThreads(threads = 16)
addArchRGenome("hg38")
library(pheatmap)
`%notin%` = Negate(`%in%`)

start_time <- Sys.time()

# Create directory for plots
plots_dir = "analysis_plots/"
dir.create(plots_dir)

current_dir <- getwd()

# load in the project from step_1
working_dir <- paste0(getwd(), "/ArchR_project/")
proj2 <- loadArchRProject(path = working_dir,
                          force = FALSE,
                          showLogo = FALSE)


# produce plot of TSS and nFrags prior to enrichment
df <- getCellColData(proj2, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(x = df[,1], y = df[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment", xlim = c(log10(500), quantile(df[,1], probs = 0.99)), ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 3, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
ggsave(p, filename = paste0(plots_dir, "pre_filtering_TSS_v_nFrag.pdf")) # save the plot

# filter out to highly quality nuclei
idxPass <- which(proj2$nFrags >= 10^3.2 & proj2$TSSEnrichment >= 5)
cellsPass <- proj2$cellNames[idxPass]
proj2 <- proj2[cellsPass, ]

df <- getCellColData(proj2, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(x = df[,1], y = df[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment", xlim = c(log10(500), quantile(df[,1], probs = 0.99)), ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 3, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
ggsave(p, filename = paste0(plots_dir, "post_filtering_TSS_v_nFrag.pdf")) # save the plot

# Now, actually filter doublets from ArchR project
proj2 <- filterDoublets(proj2, filterRatio = 2)

print("Perform Leiden clustering")
# Perform iterative LSI dimensionality reduction
proj2 <- addIterativeLSI(ArchRProj = proj2, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE,
                         iterations = 5, clusterParams = list(resolution = 0.2, sampleCells = 10000, maxClusters = 15),
                         varFeatures = 10000)


print("Perform batch correction")
# Perform batch effect correction with Harmony
proj2 <- addHarmony(ArchRProj = proj2, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample", force = TRUE)

print("Perform UMAP")
# Perform UMAP
proj2 <- addUMAP(ArchRProj = proj2, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors = 30,
                 minDist = 0.5, metric = "cosine", force = TRUE)

# Perform clustering -- do not do this as it takes a very long time
#proj2 <- addClusters(input = proj2, reducedDims = "Harmony", method = "Seurat", name = "HarmonyClusters", resolution = 1, force = TRUE,
#                    filterBias=TRUE)

# save the harmonized project 
saveArchRProject(ArchRProj = proj2,
                 outputDirectory = paste0(current_dir, "/harmonized_ArchRproject"), load = FALSE)

print("Analysis step 2 done!")

end_time <- Sys.time()
time_taken <- as.numeric(difftime(end_time, start_time, units = "hours"))
# Print the total time taken 
cat("Total time taken:", time_taken, "hours\n")
