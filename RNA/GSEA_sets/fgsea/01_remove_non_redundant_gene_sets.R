# Many of the gene sets in gene set enrichment analysis are very similar
# To improve interpretability, we will try to simplify the gene sets into non-redundant ones 
# There is a simplify() function in clusterProfiler, but it apparently no longer works
# Therefore, we will create a script that examines the Jaccard similarity between each of the gene sets in the Biology Process gene set and return non-redundant gene sets 

library(fgsea)
library(tidyverse)
library(pheatmap)

start_time = proc.time()

# extract the pathways
pathways <- gmtPathways("c5.go.bp.v2023.2.Hs.symbols.gmt")

# function to compute Jaccard similarity, which the length of the intersection divided by the union
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# create a similarity matrix (n x n for n pathways) to store the Jaccard similarities
pathway_names <- names(pathways)
n <- length(pathway_names)
similarity_matrix <- matrix(0, n, n, dimnames = list(pathway_names, pathway_names))

print(paste0("number of pathways in gene set", n))

# fill the matrix with Jaccard similarities
for (i in 1:n) {
    if (i %% 1000 == 0) {
        print(i)
        flush.console()
        }
    for (j in i:n) {  # iterate only for j > i
    sim <- jaccard_similarity(pathways[[i]], pathways[[j]])
    similarity_matrix[i, j] <- sim
    similarity_matrix[j, i] <- sim  # distance matrix is symmetric
  }
}

# convert similarity matrix to a distance matrix
distance_matrix <- 1 - similarity_matrix  # Jaccard similarity to distance

# perform hierarchical clustering
hc <- hclust(as.dist(distance_matrix), method = "complete")

# save these as objects 
saveRDS(distance_matrix, file = "fgsea_gene_set_distance_matrix.rds")
saveRDS(hc, file = "fgsea_hierarchical_clustering.rds")

end_time = proc.time()
elapsed_time = end_time - start_time

print(paste0("Elapsed time for this script in seconds is: ", elapsed_time))
print("Script complete!")
