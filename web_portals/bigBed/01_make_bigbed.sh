#!/bin/bash

# convert bed (fragment file) to bigBed

usage() {
    echo "Usage: $0 [-b bed_file]"
    echo "  -b  Specify the input bed file"
    exit 1
}

####
# Parse command-line options
while getopts ":b:" opt; do
    case ${opt} in
        b)
            bed_file=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

if [[ -z $bed_file ]]; then
	usage
fi

start_time=$(date +%s)

hg38_chrom_sizes="/mnt/data1/william/human_heart_project/Final_manuscript_analysis/human_genome/SCENICplus/hg38.chrom.sizes"

bed_file_name=$(basename $bed_file .bed.gz)
bed_file_copy="${bed_file_name}.bed.gz"

cp -r $bed_file $PWD/${bed_file_copy}
zcat $bed_file_copy | sort -k1,1 -k2,2n  > ${bed_file_name}.bed
bedToBigBed ${bed_file_name}.bed $hg38_chrom_sizes ${bed_file_name}.bb

# remove the intermediate files
rm $bed_file_copy
rm ${bed_file_name}.bed

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

echo "Elapsed time: ${hours} hours, ${minutes} minutes, ${seconds} seconds"
