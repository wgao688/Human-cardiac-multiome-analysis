#!/bin/bash

OUT_DIR="${PWD}"
CBDIR="${OUT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="${OUT_DIR}/human_cardiac_peaks.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"
DATABASE_PREFIX="human_cardiac"
SCRIPT_DIR="${PWD}"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 32
