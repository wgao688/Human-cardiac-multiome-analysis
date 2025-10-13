#!/bin/bash

# directory w/ 1KG VCFs
VCF_DIR="1KG_vcfs"
EUR_SAMPLES="1K_genomes_eur_samples.txt"

OUT_DIR="1KG_vcfs_EUR"
mkdir -p $OUT_DIR

vcf_files=($VCF_DIR/*.vcf.gz)
batch_size=4
total=${#vcf_files[@]}

for ((i=0; i<$total; i+=batch_size)); do
    echo "Processing batch $((i/batch_size + 1))..."
    # Start background jobs for this batch
    for ((j=i; j<i+batch_size && j<total; j++)); do
        vcf=${vcf_files[j]}
        filename=$(basename $vcf)
        out_vcf="$OUT_DIR/eur_${filename}"
        echo "Filtering $vcf for European samples..."
        bcftools view -S $EUR_SAMPLES --force-samples -Oz -o $out_vcf $vcf &
    done
    # Wait for all 4 background jobs to finish
    wait
    # Index the output files
    for ((j=i; j<i+batch_size && j<total; j++)); do
        vcf=${vcf_files[j]}
        filename=$(basename $vcf)
        out_vcf="$OUT_DIR/eur_${filename}"
        bcftools index $out_vcf
    done
done
