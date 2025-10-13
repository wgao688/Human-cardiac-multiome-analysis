#!/bin/bash

REGION_BED="../consensus_peaks.bed"

HG_DIR="/mnt/data1/william/human_heart_project/Final_manuscript_analysis/human_genome/SCENICplus/"
GENOME_FASTA="$HG_DIR/hg38.fa"
CHROMSIZES="$HG_DIR/hg38.chrom.sizes"
DATABASE_PREFIX="human_cardiac"

bash create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        human_cardiac_peaks.fa \
        1000 \
        yes
