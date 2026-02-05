import sys
sys.path.append("trajectory_funcs.py") 
import trajectory_funcs as funcs
import scanpy as sc
import scipy
from scipy.cluster.hierarchy import fcluster
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from scipy import sparse
from scipy.stats import spearmanr
import statsmodels.api as sm
import gseapy as gp
from gseapy import enrichr
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
import os
import argparse
import time
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

start_time = time.time()

parser = argparse.ArgumentParser(description="Process trajectory data for specific cell type.")
parser.add_argument("-c", "--cell_type", type=str,  required=True, help="cell_type")
args = parser.parse_args()

cell_type = cell_type=args.cell_type

main_output_dir = "trajectory_res/"
os.makedirs(main_output_dir, exist_ok=True)
output_dir = main_output_dir + "/" + cell_type + "/"
os.makedirs(output_dir, exist_ok=True)

cell_adata_dir = "../subclustering_analysis/post_subcluster_annotation/"

if cell_type == "Fibroblast" or cell_type == "Myeloid":
    adata = sc.read_h5ad(cell_adata_dir + cell_type + ".h5ad")
else:
    adata = sc.read_h5ad(cell_adata_dir + cell_type + ".h5ad")

adata_metadata = pd.read_csv(cell_adata_dir + cell_type +  "_annotation.csv", index_col=0)
adata.obs = adata_metadata

# compare non-diseased to these groups
other_groups = ['fetal:ND', 'postnatal:DCM', 'postnatal:HCM', 'postnatal:ICM']

for other_group in other_groups:
    print(other_groups, flush=True)

    # reverse order for fetal --> ND
    if other_group == "fetal:ND":
        df_results, ora_res_df, fig = funcs.run_trajectory_pipeline(adata=adata, other_group="postnatal:ND", non_diseased_group = other_group)
    else:
        df_results, ora_res_df, fig = funcs.run_trajectory_pipeline(adata=adata, other_group=other_group)

    # save the results
    df_results.to_csv(output_dir + cell_type + "_" + other_group + "_gene_correlations.csv")
    plt.show()
    plt.savefig(output_dir + cell_type + "_" + other_group + "_trajectory_heatmap.pdf", dpi=300)
    #plt.savefig(output_dir + cell_type + "_" + other_group + "_trajectory_heatmap.svg", dpi=300)
    ora_res_df.to_csv(output_dir + cell_type + "_" + other_group + "_ora_res.csv")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time is {elapsed_time} s", flush=True)
