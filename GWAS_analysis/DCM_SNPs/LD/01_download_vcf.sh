#!/bin/bash

# base URL for high-coverage 1000 Genomes VCF files
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"

vcf_dir="1KG_vcfs/"
mkdir -p $vcf_dir

cd $vcf_dir

# loop through chromosomes 1 to 22
for CHR in {1..22}; do
  
	VCF_FILE="1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  	FILE_URL="${BASE_URL}/${VCF_FILE}"

  	echo "Downloading ${VCF_FILE}..."
  	wget -q --show-progress ${FILE_URL}

  	# Full URL for the index file (.tbi)
	TBI_FILE="${VCF_FILE}.tbi"
	TBI_URL="${BASE_URL}/${TBI_FILE}"
	echo "Downloading ${TBI_FILE}..."
	wget -q --show-progress ${TBI_URL}
done

echo "All downloads completed!"
