#!/bin/bash

# Create null snATAC-seq peak sets, with n such repetitions

n=1000
# the actual peaks called by MACS3
input_bed="01_snATAC_enhancer_peaks.bed"
# the human genome file for bedtools shuffle
genome_file="human.hg38.genome"
# the GWAS SNPs
enhancer_file="Spurrell_healthy_enhancers.bed"
# the output file with the number of GWAS SNPs intersecting the null peaks
output_file="null_dist_results.txt"

# ensure the output file is empty before writing
> "$output_file"

start_time=$(date +%s)

# loop for n iterations
for i in $(seq 1 $n); do

	echo "$i"

   	# shuffle peaks (a null set) and save to a temporary file
   	bedtools shuffle -i "$input_bed" -g "$genome_file" > shuffled_snATAC_peak.bed
	
	# get the intersection of the enhancer peaks and the shuffled snATAC peak
	bedtools intersect -a $enhancer_file -b shuffled_snATAC_peak.bed -c > shuffled_enhancer_snATAC_peak_intersection.bed

    	# intersect with GWAS loci and count the overlaps
	overlap_count=$(awk '$NF !=0 {print}' shuffled_enhancer_snATAC_peak_intersection.bed | wc -l)

    	# append the result to the output file
    	echo "$overlap_count" >> "$output_file"
done

# removed temporary files
rm shuffled_snATAC_peak.bed

end_time=$(date +%s)
elapsed_time=$((end_time-start_time))

# Display elapsed time
echo "Script completed in $elapsed_time seconds."
