# Preprocesses all of the fragment files and saves to an initial ArchR project directory

library(ArchR)
library("unixtools")
set.tempdir("/mnt/data1/william/tmp")
addArchRThreads(threads = 16) 
addArchRGenome("hg38")

start_time <- Sys.time()


current_dir <- getwd()
current_dir

sample_metadata <- read.csv("01_ArchR_sample_metadata.csv", header=TRUE)
fragment_files = paste0("../donor_fragment_files/", sample_metadata$fragment_file)

print("Creating arrow files...")
# create ArrowFiles
ArrowFiles <- createArrowFiles(
  inputFiles = fragment_files,
  sampleNames = sample_metadata$sample_id,
  minTSS = 4, # Don't set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE)

print("Finishing creating arrow files")

# compute doublet scores but do not filter yet
doubScores <- addDoubletScores(input = ArrowFiles, k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1)

# create ArchR project
# use copyArrows=TRUE so that you maintain an unaltered copy for later usage.
proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = paste0(current_dir, "/ArchR_project"), copyArrows = TRUE)

# save the ArchR project
saveArchRProject(ArchRProj = proj, outputDirectory = paste0(current_dir, "/ArchR_project"), load = FALSE)

print("Step 1 analysis is complete!")

end_time <- Sys.time()
time_taken <- as.numeric(difftime(end_time, start_time, units = "hours"))
# Print the total time taken
cat("Total time taken:", time_taken, "hours\n")
