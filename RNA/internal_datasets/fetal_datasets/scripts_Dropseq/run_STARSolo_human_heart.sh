#!/bin/bash

# Function to display usage information
usage() {
  echo "Usage: $0 [-g <genomeDir>] [-t <threads>] [-1 <fastq_R1>] [-2 <fastq_R2>] [-v <technology>] [-h]"
  echo "  -g 	Path to the STAR genome index directory"
  echo "  -t 	Number of threads (default: 8)"
  echo "  -1 	Path to R1 FASTQ file (cell barcodes and UMIs)"
  echo "  -2 	Path to R2 FASTQ file (RNA reads)"
  echo "  -v 	Specify the technology (must be one of Dropseq, 10X_3p_v2, 10X_3p_v3, 10X_5p_v1, or 10X_Multiome"
  echo "  -o 	Output directory"
  echo "  -h 	Show this help message and exit"
  exit 1
}

# Default values
THREADS=8

# Parse command-line options
while getopts ":g:t:1:2:v:o:h" opt; do
  case "${opt}" in
    g)
      GENOME_DIR=${OPTARG}
      ;;
    o)
      OUTPUT_DIR=${OPTARG}
      ;;
    t)
      THREADS=${OPTARG}
      ;;
    1)
      FASTQ_R1=${OPTARG}
      ;;
    2)
      FASTQ_R2=${OPTARG}
      ;;
    v) 
      TECHNOLOGY=${OPTARG}
      ;;
    h)
      usage
      ;;
    *)
      usage
      ;;
  esac
done

# Check that mandatory arguments are provided
if [[ -z "$GENOME_DIR" || -z "$FASTQ_R1" || -z "$FASTQ_R2" || -z "$TECHNOLOGY" ]]; then
    echo "Error: Missing required argument(s)."
    usage
fi

echo "Arguments received:"
echo "GENOME_DIR: $GENOME_DIR"
echo "THREADS: $THREADS"
echo "FASTQ_R1: $FASTQ_R1"
echo "FASTQ_R2: $FASTQ_R2"
echo "TECHNOLOGY: $TECHNOLOGY"

OUTPUT_DIR=$(dirname $FASTQ_R1)/$OUTPUT_DIR
echo "OUTPUT_DIR: $OUTPUT_DIR"

#rm -rf $OUTPUT_DIR

# adjust the parameters for --soloUMIlen, --soloCBlen and --soloCBwhitelist according to the technology
# for Dropseq: --soloCBlen 12 --soloUMIlen 8 --soloCBwhitelist None
# for 10X_3p_v2: --soloCBlen 16 --soloUMIlen 10 --soloCBwhitelist: 2016 10x barcodes: /home/william/cellranger-8.0.0/lib/python/cellranger/barcodes/737K-august-2016.txt  
# for 10X_3p_v1: --soloCBlen 16 --soloUMIlen 10 --soloCBwhitelist: 2018 10x barcodes: /home/william/cellranger-8.0.0/lib/python/cellranger/barcodes/3M-february-2018.txt   
# for 10X_5p_v1: --soloCBlen 16 --soloUMIlen 10 --soloCBwhitelist: 2016 10x barcodes: /home/william/cellranger-8.0.0/lib/python/cellranger/barcodes/737K-august-2016.txt  
# for 10X_multiome: --soloCBlen 16 --soloUMI 12 --soloCBwhitelist: Multiome 10x barcodes: /home/william/cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt

# the technology imputted must be one of the following
if [[ "$TECHNOLOGY" == "Dropseq" ]]; then
	CB_LEN=12
	UMI_LEN=8
	WHITELIST="None"
elif [[ "$TECHNOLOGY" == "10X_3p_v2" ]]; then
	CB_LEN=16
        UMI_LEN=10
        WHITELIST="/home/william/cellranger-8.0.0/lib/python/cellranger/barcodes/737K-august-2016.txt"
elif [[ "$TECHNOLOGY" == "10X_3p_v3" ]]; then
	CB_LEN=16
        UMI_LEN=12
        WHITELIST="/home/william/cellranger-8.0.0/lib/python/cellranger/barcodes/3M-february-2018.txt"
elif [[ "$TECHNOLOGY" == "10X_5p_v1" ]]; then
	CB_LEN=16
        UMI_LEN=10
        WHITELIST="/home/william/cellranger-8.0.0/lib/python/cellranger/barcodes/737K-august-2016.txt"
elif [[ "$TECHNOLOGY" == "10X_Multiome" ]]; then
	CB_LEN=16
        UMI_LEN=12
        WHITELIST="/home/william/cellranger-8.0.0/lib/python/cellranger/barcodes/737K-arc-v1.txt"
else
	echo "TECHNOLOGY is not valid. Please see -h for the valid arguments."
    	exit 1
fi

# Run STARsolo, load read2 first (sequence) and then read1 UMI)
STAR --runThreadN $THREADS \
     --genomeDir $GENOME_DIR \
     --readFilesIn $FASTQ_R2 $FASTQ_R1 \
     --readFilesCommand zcat \
     --soloType CB_UMI_Simple \
     --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
     --soloCBwhitelist $WHITELIST \
     --soloUMIlen $UMI_LEN \
     --soloCBlen $CB_LEN \
     --soloFeatures Gene GeneFull \
     --outFileNamePrefix $OUTPUT_DIR/ \
     --outTmpDir $OUTPUT_DIR/tmp/ \
     --soloCellFilter EmptyDrops_CR \
     --limitBAMsortRAM 70000000000 --outSAMtype BAM SortedByCoordinate

rm $OUTPUT_DIR/Aligned.out.sam
echo "STARSolo alignment and quantification complete."
