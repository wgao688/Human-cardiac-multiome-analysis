# Human cardiac multiome analysis
This repository contains the scripts used for the analysis of the human cardiac snRNA-seq, snATAC-seq, and spatial transcriptomic datasets in the manuscript titled "Single-cell multiomic integration identifies widespread, cell-type resolved fetal reactivation in the diseased human heart"). For any questions regarding these scripts, please reach out to william.gao@pennmedicine.upenn.edu. 

This repository includes several subdirectories. For each subdirectory, there is a `00_info.txt` file describing the contents of the directory, and the order of the commands used to generate the paper figures. Some of the scripts are run in the command line; others are interactive Jupyter notebooks. 

To run these notebooks/scripts, there are some R/python libraries that need to be downloaded, which will be indicated in the scripts/notebooks. Because some files are too large, they are not included in this directory, but can be reproduced using the scripts, as described in the `00_info.txt` files within each directory. These large files include *.h5ad (anndata format) objects with snRNA-seq + snATAC-seq + spatial objects, and FASTQ files. For example, for any directories that specify raw data, please download these either from the external links for datasets from other studies, or GEO for the datasets generated in this study. The new samples generated as part of this study (19 snRNA-seq donors, 11 snATAC-seq donors) are deposited as raw FASTQ files in the Gene Expression Omnibus under the accession number *TBD*. 

Here is the structure of this directory. The snATAC-seq directory is `ATAC`, the snRNA-seq and spatial RNA-seq directory is `RNA/`. 
```
├── 00_info.txt
├── ATAC
│   ├── aggregated_analysis
│   │   ├── 00_info.txt
│   │   ├── 00_original_metadata.csv
│   │   ├── 00_update_metadata.ipynb
│   │   ├── 00_updated_metadata.csv
│   │   ├── 01A_examine_QC.py
│   │   ├── 01B_QC_filter_per_sample.txt
│   │   ├── 01C_QC_updated_metadata.csv
│   │   ├── 01C_update_metadata_with_QC.ipynb
│   │   ├── 01D_filter_by_QC.py
│   │   ├── 02_combine_filtered_adatas.py
│   │   ├── 03_harmonize_adata.py
│   │   ├── 04B_annotate_leiden_clusters.ipynb
│   │   ├── 04_impute_gene_expression.py
│   │   ├── 05A_clean_up_peak_adata.py
│   │   ├── 05B_subset_to_individual_cell_types.py
│   │   ├── 05C_annotate_indiv_cell_type.ipynb
│   │   ├── 05D_reannotate_cell_types.py
│   │   ├── 06B_merge_peaks.py
│   │   ├── 06C_run_peak_UMAP.py
│   │   ├── 06D_examine_UMAP_and_obtain_marker_peaks.ipynb
│   │   ├── 06E_SnapATAC2_motifs_and_marker_peaks.ipynb
│   │   ├── 06_call_peaks.py
│   │   ├── 07A_make_ATAC_adata_consistent_with_RNA.ipynb
│   │   ├── 07B_get_peak_bed.ipynb
│   │   ├── 07B_snATAC_peaks.bed
│   │   └── 07C_annotate_peaks_ChipSeeker.ipynb
│   ├── external_datasets
│   │   ├── 00_info.txt
│   │   ├── Ameen
│   │   ├── ENCODE
│   │   ├── Kanemaru
│   │   └── Kuppe
│   └── internal_datasets
│       ├── 00_directories.txt
│       ├── 00_info.txt
│       ├── 00_send_cellranger_atac.sh
│       ├── 01_get_fragment_files.sh
│       ├── 01_metadata.txt
│       ├── 01_visualize_metadata.ipynb
│       ├── 02_full_metadata.txt
│       ├── 02_verify_full_metadata.ipynb
│       ├── 03_send_cellsnplite.sh
│       ├── 04_verify_sex_BAM.ipynb
│       └── run_cellsnplite_pseudobulked.sh
├── README.md
├── RNA
│   ├── 00_info.txt
│   ├── GSEA_sets
│   │   └── fgsea
│   │       ├── 00_info.txt
│   │       ├── 02_examine_gene_set_hierarchical_clusters.ipynb
│   │       ├── Step_01_remove_non_redundant_gene_sets.R
│   │       ├── fgseapy_at_0.7_dissimilarity_GO_BP_terms.csv
│   │       └── filtered_GO_Hs_symbols.gmt
│   ├── aggregated_analysis
│   │   ├── 00_info.txt
│   │   ├── 01A_combine_LV_ND_adata.ipynb
│   │   ├── 01B_combine_fetal_LV.ipynb
│   │   ├── 01C_combine_LV_D.ipynb
│   │   ├── 02_combine_all_datasets.py
│   │   ├── 03_check_metadata_consistency.ipynb
│   │   ├── 04A_run_harmony_integration.py
│   │   ├── 04B_visualize_and_annotate_post_harmony.ipynb
│   │   ├── 04C_run_no_integration.py
│   │   ├── 05A_run_scVI.py
│   │   ├── 05C_run_subclustering.py
│   │   ├── 05E_combine_subclustered_adata.ipynb
│   │   ├── 06A_visualize_UMAP_embeddings.ipynb
│   │   ├── 06B_LISI_and_runtime_computation.ipynb
│   │   ├── 06_adata_metadata_for_LISI.csv
│   │   ├── 06_integration_runtimes.txt
│   │   ├── 07A_make_RNA_metadata_consistent.ipynb
│   │   ├── 07B_subsample_adata.ipynb
│   │   ├── cell_cell_communication
│   │   │   └── liana
│   │   │       ├── 00_info.txt
│   │   │       ├── 01_run_liana_for_heart_dataset_part_1.py
│   │   │       ├── 02_run_liana_tensor_decomposition.ipynb
│   │   │       ├── 03_analysis_post_tensor_decomposition.ipynb
│   │   │       └── 04_investigate_TGF_beta_loadings.ipynb
│   │   ├── cell_type_proportion_analysis
│   │   ├── compare_integration_methods
│   │   ├── create_sce
│   │   ├── metadata_plots
│   │   │   ├── 00_info.txt
│   │   │   ├── 01_create_donor_level_metadata.ipynb
│   │   │   ├── 01_donor_level_metadata.csv
│   │   │   ├── 02_ggplot_donor_metadata.ipynb
│   │   │   ├── 03_examine_annotation_concordance.ipynb
│   │   │   └── 04_produce_snRNA_UMAP_dotplots.ipynb
│   │   ├── original_annotation_vs_revised
│   │   │   ├── 00_info.txt
│   │   │   └── 01_compare_original_vs_revised.py
│   │   ├── senescence_analysis
│   │   │   ├── 00_extract_Senescence_gene_set_list.ipynb
│   │   │   ├── 00_info.txt
│   │   │   ├── 01_senescence_analysis.ipynb
│   │   │   ├── 01_senescence_results.csv
│   │   │   ├── 02_produce_senescence_analysis_plots.ipynb
│   │   │   ├── R-HSA-2559583_list.txt
│   │   │   ├── SenMayo_list.txt
│   │   │   └── Senescence_Gene_Set_dict.txt
│   │   └── transcriptional_variability_analysis
│   │       └── decibel
│   │           ├── 00_info.txt
│   │           ├── 01A_run_decibel_with_subsampling.py
│   │           ├── 01B_run_decibel_without_subsampling.py
│   │           ├── 01_all_decibel_results.csv
│   │           └── 02_analyze_decibel_subsampled_results.ipynb
│   ├── external_datasets
│   │   ├── Chaffin_et_al_2022
│   │   │   ├── 00_info.txt
│   │   │   ├── 01_preprocess_Chaffin_non_diseased.ipynb
│   │   │   └── 02_preprocess_Chaffin_diseased.ipynb
│   │   ├── ENCODE_v4
│   │   │   ├── 00_info.txt
│   │   │   ├── 01_examine_and_filter_metadata.ipynb
│   │   │   ├── 02_create_combined_count_matrix.ipynb
│   │   │   ├── 03_preprocess_combined_adata.ipynb
│   │   │   ├── 04_run_scrublet.py
│   │   │   ├── ENCODE_files.txt
│   │   │   ├── filtered_metadata.tsv
│   │   │   └── metadata.tsv
│   │   ├── Heart_Atlas_v2
│   │   ├── Hill_et_al_2022
│   │   │   ├── 00_info.txt
│   │   │   ├── 01_load_Hill_to_adata.ipynb
│   │   │   ├── 02A_preprocess_Hill_ND_adata.ipynb
│   │   │   └── 02B_preprocess_Hill_D_adata.ipynb
│   │   ├── Koenig_et_al_2022
│   │   │   ├── 00_donor_metadata.txt
│   │   │   ├── 00_info.txt
│   │   │   ├── 01_convert_Seurat_to_adata.R
│   │   │   ├── 02_load_into_adata.py
│   │   │   ├── 03_reformat_all_Koenig_adata.ipynb
│   │   │   ├── 04A_preprocess_Koenig_nuclei_non_diseased.ipynb
│   │   │   ├── 04B_preprocess_Koenig_nuclei_diseased.ipynb
│   │   │   └── 04C_preprocess_Koenig_cell.ipynb
│   │   ├── Kuppe_et_al_2022
│   │   │   ├── 00_healthy_donors.txt
│   │   │   ├── 00_info.txt
│   │   │   ├── 01_inspect_all_Kuppe_adata.ipynb
│   │   │   ├── 02A_preprocess_Kuppe_ND.ipynb
│   │   │   └── 02B_preprocess_Kuppe_D.ipynb
│   │   ├── Reichart_et_al_2022
│   │   │   ├── 00_info.txt
│   │   │   ├── 01A_preprocess_Reichart_LV_ND.ipynb
│   │   │   ├── 02_preprocess_Reichart_LV_diseased.ipynb
│   │   │   └── 02_preprocess_Reichart_LV_diseased.py
│   │   ├── Simonson_et_al_2023
│   │   │   ├── 00_info.txt
│   │   │   ├── preprocess_Simonson_diseased.ipynb
│   │   │   └── preprocess_Simonson_non_diseased.ipynb
│   │   └── fetal_datasets
│   │       ├── ENCODE_v4
│   │       │   ├── 00_download_tarballs.sh
│   │       │   ├── 00_fetal_multiome_metadata.txt
│   │       │   ├── 00_info.txt
│   │       │   ├── 01A_unpack_tarballs.sh
│   │       │   ├── 01B_gunzip_files.sh
│   │       │   ├── 02_generate_combined_adata.ipynb
│   │       │   └── 03_reformat_adata.ipynb
│   │       └── Sim_et_al_2020
│   │           ├── 00_process_tar.py
│   │           ├── 01_preprocess_all_Sim.ipynb
│   │           ├── 02A_preprocess_Sim_fetal.ipynb
│   │           └── 02B_preprocess_Sim_postnatal_ND.ipynb
│   ├── internal_datasets
│   │   ├── fetal_datasets
│   │   │   ├── 00_info.txt
│   │   │   ├── 01_fetal_updated_metadata.csv
│   │   │   ├── 01_view_cellsnplite_results.ipynb
│   │   │   ├── combined_analysis
│   │   │   │   ├── 01_combine_SoupX_corrected_count_matrices.ipynb
│   │   │   │   ├── 02_scanpy_combine_adata.ipynb
│   │   │   │   └── 03_preprocess_Penn_fetal.ipynb
│   │   │   ├── fetal_datasets
│   │   │   │   ├── 00_info.txt
│   │   │   │   ├── 01_fetal_updated_metadata.csv
│   │   │   │   └── 01_view_cellsnplite_results.ipynb
│   │   │   ├── raw_data_Dropseq
│   │   │   └── scripts_Dropseq
│   │   │       ├── 00_create_directories.sh
│   │   │       ├── 00_history.txt
│   │   │       ├── 01_run_STARSolo_parallel.sh
│   │   │       ├── 02_index_bam_files.sh
│   │   │       ├── 03_gzip_STAR_files.sh
│   │   │       ├── 04_send_SoupX.sh
│   │   │       ├── 05_send_cellsnplite.sh
│   │   │       ├── run_STARSolo_human_heart.sh
│   │   │       ├── run_SoupX.R
│   │   │       └── run_cellsnplite_pseudobulked.sh
│   │   └── postnatal_datasets
│   │       ├── 00_info.txt
│   │       ├── 01_metadata.txt
│   │       ├── combined_analysis
│   │       │   ├── 01_combine_SoupX_corrected_count_matrices.ipynb
│   │       │   ├── 02_scanpy_combine_adata.ipynb
│   │       │   ├── 03_preprocess_Penn_diseased_LV.ipynb
│   │       │   └── 03_preprocess_Penn_non_diseased_LV.ipynb
│   │       ├── raw_data
│   │       └── scripts
│   │           ├── 00_create_directories.sh
│   │           ├── 01_run_STARSolo_parallel.sh
│   │           ├── 02_index_bam_files.sh
│   │           ├── 03_gzip_STAR_files.sh
│   │           ├── run_STARSolo_human_heart.sh
│   │           ├── run_SoupX.R
│   │           └── send_SoupX.sh
│   └── spatial_RNA
└── human_genome
    ├── 00_info.txt
    ├── hg38-blacklist.v2.bed
    ├── hg38.blacklist.bed
    └── tss
        ├── 00_info.txt
        ├── 01_reformat_to_bed.ipynb
        ├── chromosomes.txt
        ├── tss.bed
        ├── tss.txt
        └── updated_tss.bed
```
