import gseapy as gp
from gseapy import barplot, dotplot
import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import os
import time

start_time = time.time()

def run_gseapy_for_contrast(cell_type, contrast, jaccard_gene_set):

    DESeq2_results_path = "pydeseq2_results/" + cell_type + "_" + contrast + "_results.csv"

    if not os.path.isfile(DESeq2_results_path):
        print(f"File not found: {DESeq2_results_path}. Skipping GSEA for {cell_type} - {contrast}")
        return None

    de_results = pd.read_csv(DESeq2_results_path, index_col = 0)
    gene_list = de_results[["log2FoldChange"]].sort_values(by = "log2FoldChange", ascending=False)

    # Run Gene Set Enrichment Analysis (GSEA)
    gsea_results = gp.prerank(
        rnk=gene_list,  # Your ranked list
        gene_sets=jaccard_gene_set,  # Choose the gene set, or use your own .gmt file
        outdir='gsea_results',  # Output directory
        permutation_num=1000,  # Default number of permutations
        format='png',  # Output format for plots
    )

    # extract the results as df
    gsea_results_df = gsea_results.res2d

    # remove non-finite values 
    gsea_results_df = gsea_results_df[gsea_results_df["FDR q-val"].notna() & gsea_results_df["FDR q-val"].apply(lambda x: np.isfinite(x))]

    return(gsea_results_df)

gsea_results_dir = "gseapy_results_dir/"
os.makedirs(gsea_results_dir, exist_ok=True)

#jaccard_gene_set = '../../../GSEA_sets/fgsea/filtered_GO_Hs_symbols.gmt'
jaccard_gene_set = "MSigDB_Hallmark_2020"

#cell_types = ["Adipocyte", "Cardiomyocyte", "Endocardial", "Endothelial", "Epicardial", "Fibroblast", "LEC",
#              "Lymphoid", "Mast", "Myeloid", "Pericyte", "Neuronal", "vSMC"]
cell_types = ["Cardiomyocyte", "Endocardial", "Endothelial", "Epicardial", "Fibroblast", "LEC",
              "Lymphoid", "Myeloid", "Pericyte", "Neuronal", "vSMC"]

contrast_names = ["development"]
corresponding_contrast_file_names = ["age_group_fetal_vs_young"]

for cell_type in cell_types:
    for contrast_name, contrast_file_name in zip(contrast_names, corresponding_contrast_file_names):
        
        print(f"Running gseapy for {cell_type} and contrast={contrast_name}", flush=True)
        results_df = run_gseapy_for_contrast(cell_type, contrast_file_name, jaccard_gene_set)

        if results_df is not None:
            results_df.to_csv(gsea_results_dir + cell_type + "_" + contrast_file_name + "_gseapy_results.csv")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script finished in {elapsed_time}", flush=True)
