### Perform overlap of scE2G predictions with the GWAS SNPs (lead SNPs and those within LD of r2>0.2)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import time 


start_time = time.time()

output_dir = "intersecting_SNPs/"
os.makedirs(output_dir, exist_ok=True)

def parse_chr_start_end(df, col_name, prefix=""):
    chroms = []
    starts = []
    ends = []
    for region in df[col_name]:
        chrom, pos = region.split(":")
        start, end = map(int, pos.split("-"))
        chroms.append(chrom)
        starts.append(start)
        ends.append(end)
    df[prefix + "chrom"] = chroms
    df[prefix + "start"] = starts
    df[prefix + "end"] = ends
    return df


def find_intersections(likely_links, combined_GWAS_df):
    intersections = []

    # Iterate over regions
    for idx, row in likely_links.iterrows():
        r_chrom, r_start, r_end = parse_region(row['region'])

        # Subset SNPs on the same chromosome
        snps_on_chrom = combined_GWAS_df[combined_GWAS_df['candidate_SNP'].str.startswith(r_chrom + ":")]

        for snp_idx, snp_row in snps_on_chrom.iterrows():
            s_chrom, s_pos = parse_snp(snp_row['candidate_SNP'])
            if r_start <= s_pos <= r_end:
                # Combine info from both rows
                merged_row = {**row.to_dict(), **snp_row.to_dict()}
                intersections.append(merged_row)

    return pd.DataFrame(intersections)

# cell type for prediction
cell_types = ["Cardiomyocyte", "Fibroblast", "Endothelial", "Epicardial", "Lymphoid", "Myeloid", "Neuronal", "Pericyte", "vSMC"]

# load in the GWAS results
DCM_SNPs = pd.read_csv("../GWAS/GWAS_DCM/LD/08_lead_and_LD_SNPs.csv")
DCM_SNPs['disease'] = "DCM"

filt_DCM_SNPs = DCM_SNPs[["CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B",
                          "R2", "lead_SNP", "candidate_gene", "disease", "distance", "locus"]].rename(columns = {"candidate_gene": "prioritized_gene"})

HCM_SNPs = pd.read_csv("../GWAS/GWAS_HCM/LD/04_lead_and_LD_SNPs.csv")
HCM_SNPs['locus'] = HCM_SNPs['lead_SNP']
HCM_SNPs['disease'] = "HCM"

combined_GWAS_df = pd.concat([HCM_SNPs, filt_DCM_SNPs])

# also add the candidate SNP
combined_GWAS_df['candidate_SNP'] = ("chr" + combined_GWAS_df['CHR_B'].astype(str) + ":"  +
                               combined_GWAS_df['BP_B'].astype(str) + "-" +
                               combined_GWAS_df['BP_B'].astype(str) )

# add chrom info for intersection function
combined_GWAS_df = parse_chr_start_end(combined_GWAS_df, "candidate_SNP")

threshold=0.17

for cell_type in cell_types:

    print(cell_type, flush=True)

    # load in all of the results
    e2G_df = pd.read_csv("../scE2G/results/all_multiome_only/" + cell_type + "/multiome_powerlaw_v3/encode_e2g_predictions.tsv.gz",
                       compression='gzip', delimiter="\t")

    # filter to high confidence E->P linkages
    likely_links = e2G_df[e2G_df['E2G.Score'] > threshold].copy()

    likely_links['region'] = (likely_links['chr'].astype(str) + ":"  +
                               likely_links['start'].astype(str) + "-" +
                               likely_links['end'].astype(str) )

    likely_links = parse_chr_start_end(likely_links, "region")
    combined_GWAS_df = parse_chr_start_end(combined_GWAS_df, "candidate_SNP")
    
    # find intersections
    intersections = []
    for chrom, group in likely_links.groupby("chrom"):
        snps_on_chrom = combined_GWAS_df[combined_GWAS_df["chrom"] == chrom]
        for idx, row in group.iterrows():
            overlaps = snps_on_chrom[
                (snps_on_chrom["start"] <= row["end"]) &
                (snps_on_chrom["end"] >= row["start"])
            ]
            
            if not overlaps.empty:
                for _, snp_row in overlaps.iterrows():
                    merged = {**row.to_dict(), **snp_row.to_dict()}
                    intersections.append(merged)

    # merge results together
    intersect_df = pd.DataFrame(intersections)

    intersect_df.to_csv(output_dir + cell_type + "_loci_intersections.csv", index=False)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time for script is {elapsed_time} s", flush=True)
