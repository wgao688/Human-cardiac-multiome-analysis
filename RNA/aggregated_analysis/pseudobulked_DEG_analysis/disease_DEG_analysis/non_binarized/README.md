### For disease analysis, perform non-binarized disease analysis for HCM, DCM, ICM

### STEP 1: Perform non-binarized DEG analysis.

Only perform it for cell types where there are enough donors available across multiple studies (to avoid singular matrix) for each disease subtype. Therefore, we will only examine 9 cell types: Adipocyte, Cardiomyocyte, Endothelial, Fibroblast, LEC, Lymphoid, Myeloid, Neuronal, Pericyte

Perform using this script `01_run_pydeseq2_disease_non_binary.py`

### STEP 2: Examine DEG overlaps between disease subtypes

Using these scripts `02A_examine_disease_overlap.py`, `02B_number_of_shared_DEGs.ipynb`. Perform GSEA using `02C_run_prerank_GSEA.ipynb` and visualize using `02D_visualize_GSEA_prerank.ipynb`.

### STEP 3: Examine the Jaccard similarity between DEG sets.

Using these scripts `03A_non_binarized_disease_Jaccard.ipynb` and `03B_visualize_Jaccard_similarity.ipynb`.

### STEP 4: Examine intersection between disease subtype DEGs and fetal DEGs

Using these scripts `04A_non_binarized_disease_fetal_intersect.ipynb`, `04B_produce_fetal_intersection_plot.ipynb`, and `04C_identify_percentage_of_fetalization_genes.ipynb`
