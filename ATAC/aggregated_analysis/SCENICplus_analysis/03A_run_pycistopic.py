# Run the entire pycistopic workflow
# use the scenicplus conda env

import warnings
warnings.simplefilter(action='ignore')
import pycisTopic
import pickle
import os
import time
import requests

from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import run_cgs_models_mallet
from pycisTopic.lda_models import evaluate_models

import pandas as pd

import scipy.io
from scipy.io import mmread
from scipy.sparse import csr_matrix

import pybiomart as pbm
import pyranges as pr

start_time = time.time()

print("Loading in count matrices and peak information...", flush = True)
# load in count matrix and make it a sparse matrix
count_matrix = mmread("peak_counts.mtx")
count_matrix = csr_matrix(count_matrix)

# path to blacklist
path_to_blacklist='/mnt/data1/william/human_heart_project/Final_manuscript_analysis/human_genome/hg38-blacklist.v2.bed'

# load in the cell metadata as well
cell_data = pd.read_csv("cell_metadata.csv", index_col = 0)
cell_data.index = cell_data.index.astype(str)
peak_data = pd.read_csv("peak_names.csv", index_col = 0)

# create cistopic object
print("Creating pycistopic object...", flush = True)
cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix, path_to_blacklist=path_to_blacklist, 
                                      cell_names = cell_data.index, region_names = peak_data.index)

# add cell data to a cistopic object
cistopic_obj.add_cell_data(cell_data)

print(cistopic_obj, flush = True)

out_dir = "pycistopic_outs/"
os.makedirs(out_dir, exist_ok = True)

# specify increased RAM for Mallet
os.environ['MALLET_MEMORY'] = '300G'

# Run models
mallet_path="/mnt/data1/william/SCENICplus_analysis/Mallet-202108/bin/mallet"

print("Running topic modeling with Mallet ...", flush = True)
# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[5, 10, 20, 30, 40, 50],
    n_cpu=24,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path="/mnt/data1/william/tmp/",
    save_path=out_dir,
    mallet_path=mallet_path,
)

# create a directory to store the models
os.makedirs(out_dir + 'models', exist_ok=True)

model=evaluate_models(models,
                     select_model=None, 
                     return_model=True, 
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False,
                     save = os.path.join(out_dir + 'models/model_selection.pdf'))

# add model to cisTopicObject
cistopic_obj.add_LDA_model(model)

pickle.dump(models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)

# save pycistopic object
print("Saving the models and cistopic object...", flush=True)
with open(out_dir + 'cisTopicObject.pkl', 'wb') as f:
  pickle.dump(cistopic_obj, f)

####  Run steps after topic modeling #### 

# load back in the cistopic object
cto_infile = open(out_dir + 'cisTopicObject.pkl', 'rb')
cistopic_obj = pickle.load(cto_infile)
cto_infile.close()

# load back in the models
model_infile = open(out_dir + 'models.pkl', 'rb')
models = pickle.load(model_infile)
model_infile.close()

from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)

# find clusters for the cistopic object
find_clusters(
    cistopic_obj,
    target  = 'cell',
    k = 10,
    res = [0.6, 1.2, 3],
    prefix = 'pycisTopic_',
    scale = True,
    split_pattern = '-'
)

# run UMAP
run_umap(cistopic_obj,  target  = 'cell', scale=True)

# create plots for visualizations
os.makedirs( os.path.join(out_dir,'visualization') , exist_ok = True)

# plot the UMAP
plot_metadata(cistopic_obj,
                 reduction_name='UMAP',
                 variables=['cell_type', 'pycisTopic_leiden_10_0.6'], # Labels from RNA and new clusters
                 target='cell', num_columns=3,
                 text_size=10,
                 dot_size=5,
                 figsize=(15,5),
                 save= os.path.join(out_dir, 'visualization/dimensionality_reduction_cell_type_and_leiden.pdf'))

plot_topic(cistopic_obj,
            reduction_name = 'UMAP',
            target = 'cell',
            num_columns=5,
            save = os.path.join(out_dir, 'visualization/dimensionality_reduction_topic_contr.pdf'))

# topic binarization
print("Performing binarization of topics...", flush=True)
os.makedirs(out_dir + 'topic_binarization', exist_ok = True)

# create directories to the region sets
os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok = True)

from pycisTopic.topic_binarization import *

region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=True, num_columns=5, 
    save=out_dir + '/topic_binarization/otsu.pdf'
)

region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)

binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)

from pycisTopic.utils import region_names_to_coordinates

for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

from pycisTopic.topic_qc import *
from pycisTopic.utils import fig2img

# compute topic metrics
topic_qc_metrics = compute_topic_metrics(cistopic_obj)

fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)

# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1

plt.subplots_adjust(wspace=0, hspace=-0.70)
fig.savefig(out_dir + '/topic_binarization/Topic_qc.pdf', bbox_inches='tight')

topic_annot = topic_annotation(cistopic_obj, annot_var='cell_type',
                               binarized_cell_topic=binarized_cell_topic, general_topic_thr = 0.2)

topic_qc_metrics = pd.concat([topic_annot[['cell_type', 'Ratio_cells_in_topic', 'Ratio_group_in_population']], topic_qc_metrics], axis=1)
topic_qc_metrics

# Save each of the files
with open(out_dir + '/topic_binarization/Topic_qc_metrics_annot.pkl', 'wb') as f:
  pickle.dump(topic_qc_metrics, f)
with open(out_dir + '/topic_binarization/binarized_cell_topic.pkl', 'wb') as f:
  pickle.dump(binarized_cell_topic, f)
with open(out_dir + '/topic_binarization/binarized_topic_region_otsu.pkl', 'wb') as f:
  pickle.dump(region_bin_topics_otsu, f)
with open(out_dir + '/topic_binarization/binarized_topic_region_top_3K.pkl', 'wb') as f:
  pickle.dump(region_bin_topics_top_3k, f)


# perform DAR analysis
print("Perform DAR analysis...", flush = True)

from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

os.makedirs(out_dir + '/DARs', exist_ok = True)

variable_regions = find_highly_variable_features(normalized_imputed_acc_obj,
                                           min_disp = 0.05,
                                           min_mean = 0.0125,
                                           max_mean = 3,
                                           max_disp = np.inf,
                                           n_bins=20,
                                           n_top_features=None,
                                           plot=True,
                                           save= out_dir + '/DARs/HVR_plot.pdf')
print("Number of variable regions: " + str(len(variable_regions)))

markers_dict = find_diff_features(cistopic_obj,
                      imputed_acc_obj,
                      variable='cell_type',
                      var_features=variable_regions,
                      contrasts=None,
                      adjpval_thr=0.05,
                      log2fc_thr=np.log2(1.5),
                      n_cpu=5,
                      _temp_dir= '/mnt/data1/william/tmp',
                      split_pattern = '-')

x = [print(x + ': '+ str(len(markers_dict[x]))) for x in markers_dict.keys()]

with open(out_dir + '/DARs/Imputed_accessibility.pkl', 'wb') as f:
  pickle.dump(imputed_acc_obj, f)

with open(out_dir + '/DARs/DARs.pkl', 'wb') as f:
  pickle.dump(markers_dict, f)


dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://nov2020.archive.ensembl.org/')
annot = dataset.query(attributes=['chromosome_name', 'start_position', 'end_position', 'strand', 'external_gene_name', 'transcription_start_site', 'transcript_biotype'])
annot['Chromosome/scaffold name'] = 'chr' + annot['Chromosome/scaffold name'].astype(str)
annot.columns=['Chromosome', 'Start', 'End', 'Strand', 'Gene','Transcription_Start_Site', 'Transcript_type']
annot = annot[annot.Transcript_type == 'protein_coding']
annot.Strand[annot.Strand == 1] = '+'
annot.Strand[annot.Strand == -1] = '-'
pr_annotation = pr.PyRanges(annot.dropna(axis = 0))

# Get chromosome sizes
target_url='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
chromsizes=pr.PyRanges(chromsizes)

# Gene score imputation
print("Perform gene activity imputation...", flush = True)
from pycisTopic.gene_activity import *
gene_act, weigths = get_gene_activity(imputed_acc_obj, # Region-cell probabilities
                pr_annotation, # Gene annotation
                chromsizes, # Chromosome size
                use_gene_boundaries=True, # Whether to use the whole search space or stop when encountering another gene
                upstream=[1000, 100000], # Search space upstream. The minimum means that even if there is a gene right next to it
                                      #these bp will be taken (1kbp here)
                downstream=[1000,100000], # Search space downstream
                distance_weight=True, # Whether to add a distance weight (an exponential function, the weight will decrease with distance)
                decay_rate=1, # Exponent for the distance exponential funciton (the higher the faster will be the decrease)
                extend_gene_body_upstream=10000, # Number of bp upstream immune to the distance weight (their value will be maximum for
                                      #this weight)
                extend_gene_body_downstream=500, # Number of bp downstream immune to the distance weight
                gene_size_weight=False, # Whether to add a weights based on the length of the gene
                gene_size_scale_factor='median', # Dividend to calculate the gene size weigth. Default is the median value of all genes
                                      #in the genome
                remove_promoters=False, # Whether to remove promoters when computing gene activity scores
                average_scores=True, # Whether to divide by the total number of region assigned to a gene when calculating the gene
                                      #activity score
                scale_factor=1, # Value to multiply for the final gene activity matrix
                extend_tss=[10,10], # Space to consider a promoter
                gini_weight = True, # Whether to add a gini index weigth. The more unique the region is, the higher this weight will be
                return_weights= True, # Whether to return the final weights
                project='Gene_activity') # Project name for the gene activity object

gene_markers_dict = find_diff_features(cistopic_obj,
                      gene_act,
                      variable='cell_type',
                      var_features=None,
                      contrasts=None,
                      adjpval_thr=0.05,
                      log2fc_thr=np.log2(1.5),
                      n_cpu=5,
                      _temp_dir= '/mnt/data1/william/tmp')

os.makedirs(out_dir + '/DAGs', exist_ok = True)

with open(out_dir + '/DAGs/DAGs.pkl', 'wb') as f:
  pickle.dump(gene_markers_dict, f)

with open(out_dir + '/DAGs/gene_activities.pkl', 'wb') as f:
  pickle.dump(gene_act, f)

end_time = time.time()
elapsed_time = end_time - start_time

print("Elapsed time in seconds is " + str(elapsed_time))
print("Script complete!")
