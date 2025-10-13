#!/bin/bash

start_time=$(date +%s)

# path to the DCM SNPs
GWAS_SNPs_path="../hg38_combined_unique_DCM_SNPs.csv"

full_SNP_name=($(awk '{print $1}' $GWAS_SNPs_path))
SNP_chr=($(awk '{print $1}' $GWAS_SNPs_path | sed 's/:/ /g' | awk '{print $1}'))

num_SNPs=${#full_SNP_name[@]}
echo "number of SNPs $num_SNPs"

plink_path="1KG_plink/"
results_dir="LD_results/"
mkdir -p $results_dir

for i in $(seq 0 $(($num_SNPs-1)) ); do
	SNP=${full_SNP_name[i]}
	chrom=${SNP_chr[i]}
	
	# get the full name of the SNP in plink
	corresponding_plink_file=$plink_path/$chrom.bim
	SNP_name_in_plink=$(awk -v snp="$SNP:" '$2 ~ snp {print $2}' $corresponding_plink_file)
	echo "$SNP_name_in_plink"

	if [[ -z $SNP_name_in_plink ]]; then
		echo "Could not find SNP $SNP in the plink file!"
	else
		echo "SNP $SNP found in the plink file"
		plink --bfile $plink_path/$chrom --ld-snp $SNP_name_in_plink --ld-window 99999 --r2 --ld-window-kb 1000 --out $results_dir/$SNP
	fi
done

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "All plink files generated!"
echo "Elapsed time: $elapsed seconds"
