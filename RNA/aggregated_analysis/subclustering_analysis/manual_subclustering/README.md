### Manual subclustering

The scripts in this directory perform manual subclustering for each of the 13 major cell types in the human heart. We use a leiden resolution of 0.25 as this tends to have the highest silhouette score and also results in fewer small clusters that are study-specific. Subclustering is guided by the top 50 marker genes per subcluster and a literature search to identify if any of these genes are associated with experimentally-validated functions in the given cell type.
