import scanpy as sc
import scvi
import time

start_time = time.time()

# import the adata
adata = sc.read_h5ad("../04_post_scrublet_adata.h5ad")

adata.layers["counts"] = adata.X.copy()  # preserve raw counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# use top 2000 genes
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
    )

print("Specifying and training the model...", flush=True)
scvi.model.SCVI.setup_anndata(
    adata, batch_key="directory",
    layer="counts",
)

model = scvi.model.SCVI(adata)
model.train()
model.save("scVI_model/", overwrite = True)
latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent

# run PCA then generate UMAP plots
print("Running PCA now...", flush=True)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)

# use scVI latent space for UMAP generation
SCVI_LATENT_KEY = "X_scVI"
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=0.3)

SCVI_CLUSTERS_KEY = "leiden_scVI"
sc.tl.leiden(adata, key_added=SCVI_CLUSTERS_KEY, resolution=0.5)
sc.pl.umap(adata, color=[SCVI_CLUSTERS_KEY], frameon=False)

# Obtain cluster-specific differentially expressed genes
print("Identifying marker genes for each cluster...", flush=True)
sc.tl.rank_genes_groups(adata, groupby=SCVI_CLUSTERS_KEY, method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata, groupby=SCVI_CLUSTERS_KEY, standard_scale="var", n_genes=3)

adata.write("scvi_adata.h5ad")

end_time = time.time()
elapsed_time = end_time - start_time

print("Elapsed time for script in seconds is: " + str(elapsed_time), flush=True)
print("Script complete!", flush=True)
