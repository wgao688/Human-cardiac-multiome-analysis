library(ArchR)
library("unixtools")
set.tempdir("/mnt/data1/william/tmp")
addArchRThreads(threads = 32)
addArchRGenome("hg38")
library(pheatmap)
library(BSgenome.Hsapiens.UCSC.hg38)
`%notin%` = Negate(`%in%`)

library(Seurat)
library(SingleCellExperiment)

start_time = Sys.time()

working_dir <- getwd()
proj3 <- loadArchRProject(path = working_dir,
                          force = FALSE,
                          showLogo = FALSE)

pathToMacs2 <- "/home/william/anaconda3/envs/macs2/bin/macs2"

# add group coverages to the ArchR project for the 'cell_type' grouping variable
proj3 <- addGroupCoverages(
  ArchRProj = proj3,
  groupBy = "cell_type"
)

proj3 <- addReproduciblePeakSet(
    ArchRProj = proj3,
    groupBy = "cell_type",
    pathToMacs2 = pathToMacs2
)

proj3 <- addReproduciblePeakSet(
    ArchRProj = proj3,
    groupBy = "cell_type",
    pathToMacs2 = pathToMacs2
)

# save the ArchR project in a new directory
saveArchRProject(ArchRProj = proj3)

end_time = Sys.time()

time_taken <- as.numeric(difftime(end_time, start_time, units = "hours"))
cat("Total time taken:", time_taken, "hours\n")
