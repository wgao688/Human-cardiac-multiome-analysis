#!/bin/bash

start_time=$(date +%s)

vcf_dir="1KG_vcfs_EUR/"
plink_dir="1KG_plink_EUR"
mkdir -p $plink_dir

cd $plink_dir

# loop through chromosomes 1 to 22
for CHR in {1..22}; do
	
	echo "$CHR"
	VCF_file="../$vcf_dir/eur_1kGP_high_coverage_Illumina.chr$CHR.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

	if [ -f $VCF_file ]; then
		echo "$VCF_file"
 		plink --vcf $VCF_file --make-bed --out $CHR
	else 
		echo "$VCF_file is not found!"
	fi
done

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "All plink files generated!"
echo "Elapsed time: $elapsed seconds"
