# Download tss.txt here:
https://epd.expasy.org/mga/hg38/gencode/gencode.html
https://epd.expasy.org/mga/hg38/gencode/gencode.v28.hg38.sga.gz

# convert .. to -
 sed -i 's/\.\./-/g' tss.txt 

# Reformat it a bed file using `01_reformat_to_bed.ipynb`

# remove the RP11 genes and they are not actually annotated genes
awk '$NF != "RP11" {print}' tss.bed > updated_tss.bed
