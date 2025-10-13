#!/bin/bash

start_time=$(date +%s)

# path to the DCM SNPs
GWAS_SNPs_path="../hg38_combined_unique_DCM_SNPs.txt"

plink_path="1KG_plink/"
results_dir="LD_results/"
mkdir -p $results_dir

tail -n +2 "$GWAS_SNPs_path" | while IFS=$'\t' read -r chr candidate_gene start end rsID total_score study length revised_study plink_names; do

    # Split comma-separated SNPs into array
    IFS=',' read -ra snp_list <<< "$plink_names"

    for SNP in "${snp_list[@]}"; do
        chrom=$(echo "$SNP" | cut -d':' -f1)
        corresponding_plink_file=$plink_path/$chrom.bim

        #SNP_name_in_plink=$(awk -v snp="$SNP" '$2 == snp {print $2}' $corresponding_plink_file)
	#echo "$SNP_name_in_plink"
	echo "$plink_names"

        #if [[ -z $SNP_name_in_plink ]]; then
        #    echo "Could not find SNP $SNP in $chrom.bim!"
        #else
        #    echo "SNP $SNP found in the plink file"
            plink --bfile $plink_path/$chrom \
                  --ld-snps $plink_names \
                  --r2 --ld-window-kb 1000 \
                  --out $results_dir/${SNP//:/_}
        #fi
    done
done < <(tail -n +2 $GWAS_SNPs_path)  # skip header line

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "All plink files generated!"
echo "Elapsed time: $elapsed seconds"
